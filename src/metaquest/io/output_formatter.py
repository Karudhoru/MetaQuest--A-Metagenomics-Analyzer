"""
MetaQuest Professional Output Formatter v5.0.0
==============================================

Enhanced formatting with better verbosity control and cleaner output.
"""

import sys
import time
import subprocess
import threading
import os
import re
from pathlib import Path
from datetime import datetime
from contextlib import contextmanager
from typing import Optional, Dict, List, Any, Tuple, Callable
import shutil

class Colors:
    """ANSI color codes for terminal output"""
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'
    
    @classmethod
    def disable(cls):
        """Disable colors for non-terminal outputs"""
        cls.HEADER = cls.BLUE = cls.CYAN = cls.GREEN = ''
        cls.YELLOW = cls.RED = cls.BOLD = cls.DIM = cls.UNDERLINE = cls.END = ''


class Spinner:
    """Animated spinner for long-running operations - Thread-safe version"""
    def __init__(self, message: str = "Processing"):
        self.message = message
        self._spinning = False  # Changed to private
        self._lock = threading.Lock()  # NEW: Thread safety
        self.thread = None
        self.frames = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏']
        self._frame_idx = 0  # Changed to private
        self._stdout = sys.stdout
    
    @property
    def spinning(self):
        """Thread-safe getter for spinning state"""
        with self._lock:
            return self._spinning
    
    @spinning.setter
    def spinning(self, value):
        """Thread-safe setter for spinning state"""
        with self._lock:
            self._spinning = value
    
    def _spin(self):
        """Spinner animation loop with thread-safe access"""
        while self.spinning:  # Now uses thread-safe property
            try:
                if self._stdout and not self._stdout.closed:
                    # Protect frame_idx from concurrent access
                    with self._lock:
                        frame = self.frames[self._frame_idx % len(self.frames)]
                        self._frame_idx += 1
                    
                    self._stdout.write(f'\r  {Colors.CYAN}[{frame}]{Colors.END} {self.message}...')
                    self._stdout.flush()
                else:
                    break
            except (ValueError, AttributeError, OSError):
                break
            time.sleep(0.1)
    
    def start(self):
        """Start the spinner"""
        if not self._stdout or self._stdout.closed:
            return
            
        self.spinning = True  # Uses thread-safe setter
        self.thread = threading.Thread(target=self._spin, daemon=True)
        try:
            self.thread.start()
        except RuntimeError:
            self.spinning = False
    
    def stop(self, final_message: Optional[str] = None):
        """Stop the spinner gracefully with thread-safe cleanup"""
        self.spinning = False  # Uses thread-safe setter
        
        # Wait for thread to finish
        if self.thread and self.thread.is_alive():
            self.thread.join(timeout=0.5)
        
        # Clear the spinner line with exclusive access
        try:
            if self._stdout and not self._stdout.closed:
                with self._lock:  # Ensure no concurrent writes
                    if final_message:
                        self._stdout.write(f'\r  {final_message}\n')
                    else:
                        self._stdout.write('\r' + ' ' * (len(self.message) + 20) + '\r')
                    self._stdout.flush()
        except (ValueError, AttributeError, OSError):
            pass

class ProgressBar:
    """Enhanced progress bar with visual bars"""
    def __init__(self, total: int, desc: str = "", unit: str = "it", width: int = 30):
        self.total = total if total > 0 else 1
        self.current = 0
        self.desc = desc
        self.unit = unit
        self.width = width
        self.start_time = time.time()
        self.last_update_time = 0

    def update(self, n: int = 1):
        """Update progress by n steps"""
        self.current = min(self.current + n, self.total)
        now = time.time()
        if now - self.last_update_time > 0.1 or self.current == self.total:
            self._display()
            self.last_update_time = now
            
    def update_to(self, value: int):
        """Update progress to an absolute value"""
        self.update(value - self.current)

    def _display(self):
        """Display the progress bar with visual elements"""
        percent = self.current / self.total
        filled = int(self.width * percent)
        bar = '█' * filled + '░' * (self.width - filled)

        elapsed = time.time() - self.start_time
        rate = self.current / elapsed if elapsed > 0 else 0
        remaining_s = (self.total - self.current) / rate if rate > 0 else 0
        remaining_str = self._format_time(remaining_s) if self.current < self.total else "done"

        status = (f"\r          {self.desc} {bar} {percent:3.0%}"
                  f" | {self.current:,}/{self.total:,} {self.unit}"
                  f" | ~{remaining_str} remaining")
        
        sys.stdout.write(status)
        sys.stdout.flush()
        if self.current >= self.total:
            print()

    def _format_time(self, seconds: float) -> str:
        if seconds < 60: return f"{seconds:.0f}s"
        if seconds < 3600: return f"{int(seconds // 60)}m {int(seconds % 60)}s"
        return f"{int(seconds // 3600)}h {int((seconds % 3600) // 60)}m"

    def __enter__(self): return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.current < self.total:
            self.update_to(self.total)


class TableFormatter:
    """Beautiful table formatting with box-drawing characters"""
    
    @staticmethod
    def format_table(title: str, sections: List[Dict[str, Any]], width: int = 73) -> str:
        """Format a beautiful table with sections"""
        lines = []
        
        # Top border with title
        lines.append(f"┌{'─' * (width - 2)}┐")
        title_padding = width - len(title) - 3
        lines.append(f"│ {Colors.BOLD}{title}{Colors.END}{' ' * title_padding}│")
        
        for i, section in enumerate(sections):
            # Section divider
            lines.append(f"├{'─' * (width - 2)}┤")
            
            # Section header (if provided)
            if 'header' in section and section['header']:
                header = section['header']
                header_padding = width - len(header) - 3
                lines.append(f"│ {Colors.CYAN}{header}{Colors.END}{' ' * header_padding}│")
                lines.append(f"├{'─' * (width - 2)}┤")
            
            # Section rows
            for key, value in section.get('rows', {}).items():
                # Skip rows with None or 0/100 values unless they're meaningful
                if value is None or (isinstance(value, str) and value == "0/100"):
                    continue
                    
                key_display = f"  {key}"
                value_display = str(value)
                
                # Calculate spacing
                key_visible = len(key_display)
                value_visible = len(value_display)
                space_needed = width - key_visible - value_visible - 3
                
                if space_needed < 1:
                    space_needed = 1
                
                lines.append(f"│{key_display}{' ' * space_needed}{value_display} │")
        
        # Bottom border
        lines.append(f"└{'─' * (width - 2)}┘")
        
        return '\n'.join(lines)
    
    @staticmethod
    def format_visual_bar(percentage: float, width: int = 10, label: str = "") -> str:
        """Create a visual percentage bar"""
        filled = int(width * (percentage / 100))
        bar = '█' * filled + '░' * (width - filled)
        
        # Color based on percentage
        if percentage >= 90:
            color = Colors.GREEN
        elif percentage >= 70:
            color = Colors.YELLOW
        else:
            color = Colors.RED
        
        result = f"{bar} {color}{percentage:5.1f}%{Colors.END}"
        if label:
            result += f" {label}"
        
        return result


class OutputFormatter:
    """Enhanced output formatter with better verbosity control"""
    
    SILENT = 0      # No output except errors
    STANDARD = 1    # Normal user output
    VERBOSE = 2     # Additional details
    DEBUG = 3       # Everything including debug

    def __init__(self, verbosity: str = 'standard', log_file: Optional[Path] = None):
        self.verbosity_map = {
            'silent': self.SILENT,
            'standard': self.STANDARD,
            'verbose': self.VERBOSE,
            'debug': self.DEBUG
        }
        self.verbosity = self.verbosity_map.get(verbosity, self.STANDARD)
        self.colors_enabled = sys.stdout.isatty()
        self.log_file = log_file
        self.start_time = time.time()
        self.stage_times = {}
        self.current_stage = None
        self._operation_stack = []
        self.is_spinning = False
        self._suppressed = False
        self._active_spinner = None  # Track current spinner for debug output pausing
        
        if not self.colors_enabled:
            Colors.disable()
        
        if self.log_file:
            self.log_file.parent.mkdir(parents=True, exist_ok=True)
            self.log_handle = open(self.log_file, 'w', encoding='utf-8')
        else:
            self.log_handle = None

    def __del__(self):
        if hasattr(self, 'log_handle') and self.log_handle:
            self.log_handle.close()

    def _log(self, message: str):
        """Write to log file if enabled"""
        if self.log_handle:
            timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            clean_msg = re.sub(r'\x1b\[[0-9;]*m', '', message)
            self.log_handle.write(f"[{timestamp}] {clean_msg}\n")
            self.log_handle.flush()

    def _print(self, message: str, min_verbosity: int = STANDARD):
        """Internal print with verbosity check"""
        # In DEBUG mode, never suppress _print so debug messages from sub-functions
        # called inside suppressed_output() context managers are still visible
        is_suppressed = self._suppressed and self.verbosity < self.DEBUG
        if not is_suppressed and self.verbosity >= min_verbosity:
            print(message)
        self._log(message)

    # ========================================================================
    # NEW: Context manager for suppressing output from external functions
    # ========================================================================
    
    @contextmanager
    def suppressed_output(self):
        """
        Context manager to suppress ALL output from external functions.
        Thread-safe version that properly handles cleanup.
        """
        old_suppressed = self._suppressed
        self._suppressed = True
        
        # Store original stdout/stderr
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        
        devnull_out = None
        devnull_err = None
        
        # Redirect to devnull only if not in debug mode
        if self.verbosity < self.DEBUG:
            try:
                devnull_out = open(os.devnull, 'w')
                devnull_err = open(os.devnull, 'w')
                sys.stdout = devnull_out
                sys.stderr = devnull_err
            except Exception:
                # If redirection fails, restore and continue
                sys.stdout = old_stdout
                sys.stderr = old_stderr
        
        try:
            yield
        finally:
            # Restore stdout/stderr first
            sys.stdout = old_stdout
            sys.stderr = old_stderr
            
            # Then close devnull files
            if devnull_out:
                try:
                    devnull_out.close()
                except Exception:
                    pass
            if devnull_err:
                try:
                    devnull_err.close()
                except Exception:
                    pass
            
            self._suppressed = old_suppressed

    # ========================================================================
    # SECTION HEADERS
    # ========================================================================
    
    def banner(self, title: str, version: str, tagline: str):
        """Print application banner"""
        if self.verbosity >= self.STANDARD:
            banner = f"""
{Colors.BOLD}{Colors.CYAN}╔{'═' * 73}╗
║  {title} v{version}{' ' * (69 - len(title) - len(version))}║
║  {tagline}{' ' * (71 - len(tagline))}║
║  MetaQuest Development Team{' ' * 45}║
╚{'═' * 73}╝{Colors.END}
"""
            print(banner)
            # Log banner summary with timestamp
            self._log(f"{title} v{version} - {tagline}")
            self._log(f"  {tagline}")

    def step_header(self, step_num: int, total_steps: int, title: str):
        """Print major pipeline step header"""
        self.current_stage = title
        self.stage_times[title] = time.time()
        if self.verbosity >= self.STANDARD:
            step_text = f"STEP {step_num}/{total_steps}: {title}"
            padding = 71 - len(step_text)
            # Log the step header as a clean one-liner with timestamp
            self._log(f"{'=' * 73}")
            self._log(f"  {step_text}")
            self._log(f"{'=' * 73}")
            # Print the decorated version to the terminal
            header = f"""
╔{'═' * 73}╗
║  {Colors.BOLD}{step_text}{Colors.END}{' ' * padding}║
╚{'═' * 73}╝"""
            print(header)

    def section_header(self, title: str):
        """Print subsection header"""
        if self.verbosity >= self.STANDARD:
            # Log clean timestamp-prefixed version to log file
            self._log(f"")
            self._log(f"  🔹 {title}")
            # Print decorated version to terminal only
            print(f"\n  {Colors.BOLD}{Colors.CYAN}🔹 {title}{Colors.END}")

    def format_step_start(self, title: str) -> str:
        """Format a step start string - used by validation engine"""
        return f"{Colors.BOLD}{Colors.CYAN}▶ {title}...{Colors.END}"

    def format_step_complete(self, title: str) -> str:
        """Format a step complete string - used by validation engine"""
        return f"{Colors.GREEN}✓ {title} validated successfully{Colors.END}"

    # ========================================================================
    # STATUS MESSAGES - Simplified tree structure
    # ========================================================================
    
    def operation(self, message: str, show_in_standard: bool = True):
        """Print operation - use sparingly, prefer spinner for long operations"""
        min_level = self.STANDARD if show_in_standard else self.VERBOSE
        if self.verbosity >= min_level:
            self._print(f"     → {message}")
        op_id = len(self._operation_stack)
        self._operation_stack.append(message)
        return op_id

    def substep(self, message: str):
        """Print sub-step with deeper tree level"""
        if self.verbosity >= self.VERBOSE:
            self._print(f"        • {message}")

    def success(self, message: str, style: Optional[str] = None):
        """Print success message"""
        if self.verbosity >= self.STANDARD:
            if style == 'bold':
                message = f"{Colors.BOLD}{message}{Colors.END}"
            self._print(f"     {Colors.GREEN}✓{Colors.END} {message}")

    def info(self, message: str, indent: int = 1):
        """Print info message"""
        if self.verbosity >= self.STANDARD:
            prefix = "     " if indent == 1 else "        "
            self._print(f"{prefix}→ {message}")

    def warning(self, message: str):
        """Print warning message"""
        if self.verbosity >= self.STANDARD:
            self._print(f"     {Colors.YELLOW}⚠{Colors.END}  {message}")

    def error(self, message: str, solutions: Optional[List[str]] = None, doc_link: Optional[str] = None):
        """Print error message with solutions"""
        if self.verbosity >= self.SILENT:  # Always show errors
            self._print(f"\n     {Colors.RED}✗ {message}{Colors.END}")
            if solutions:
                self._print(f"\n     {Colors.BOLD}💡 Solutions:{Colors.END}")
                for solution in solutions:
                    self._print(f"        • {solution}")
            if doc_link:
                self._print(f"\n     {Colors.CYAN}📚 Documentation:{Colors.END} {doc_link}")
            print()

    # ========================================================================
    # TABLE FORMATTING
    # ========================================================================
    
    def display_stats_table(self, title: str, sections: List[Dict[str, Any]], filter_empty: bool = True):
        """
        Display statistics in a beautiful table format
        
        Args:
            filter_empty: If True, skip rows with None or "0/100" values
        """
        if self.verbosity >= self.STANDARD:
            # Filter empty rows if requested
            if filter_empty:
                sections = self._filter_empty_rows(sections)
            
            table = TableFormatter.format_table(title, sections)
            self._print(f"\n{table}")
    
    def _filter_empty_rows(self, sections: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Remove empty or meaningless rows from sections"""
        filtered_sections = []
        for section in sections:
            filtered_rows = {}
            for key, value in section.get('rows', {}).items():
                # Keep non-empty meaningful values
                if value is not None and value != "0/100" and value != "N/A":
                    filtered_rows[key] = value
            
            if filtered_rows:  # Only add section if it has rows
                filtered_sections.append({
                    'header': section.get('header', ''),
                    'rows': filtered_rows
                })
        
        return filtered_sections
    
    def display_quality_metrics(self, metrics: Dict[str, float]):
        """Display quality metrics with visual bars"""
        if self.verbosity >= self.STANDARD:
            self._print("\n     📈 Quality Visualization:")
            for metric, value in metrics.items():
                metric_display = f"        {metric}"
                dots = '.' * (40 - len(metric))
                bar = TableFormatter.format_visual_bar(value, width=10, label=self._get_quality_label(value))
                self._print(f"{metric_display}{dots} {bar}")
    
    def _get_quality_label(self, percentage: float) -> str:
        """Get quality label based on percentage"""
        if percentage >= 90: return "Excellent"
        elif percentage >= 70: return "Good"
        elif percentage >= 50: return "Fair"
        else: return "Poor"

    # ========================================================================
    # PROGRESS TRACKING
    # ========================================================================
    
    @contextmanager
    def spinner(self, message: str):
        """Context manager for spinner animation"""
        if self.verbosity >= self.STANDARD and self.colors_enabled:
            spinner = Spinner(message)
            self._active_spinner = spinner
            self.is_spinning = True
            spinner.start()
            try:
                yield
            finally:
                spinner.stop()
                self.is_spinning = False
                self._active_spinner = None
        else:
            # In silent mode, just log without printing
            if self.verbosity >= self.VERBOSE:
                self.operation(message)
            yield

    def progress_bar(self, total: int, desc: str, unit: str = "it") -> ProgressBar:
        """Create a progress bar"""
        if self.verbosity >= self.STANDARD:
            return ProgressBar(total, desc, unit)
        else:
            class NoOpProgressBar:
                def update(self, n=1): pass
                def update_to(self, v): pass
                def __enter__(self): return self
                def __exit__(self, *args): pass
            return NoOpProgressBar()

    # ========================================================================
    # SUBPROCESS EXECUTION with better output control
    # ========================================================================
    
    def run_subprocess(self, cmd: List[str], operation_name: str, 
                      capture_output: bool = True, show_command: bool = False) -> Tuple[int, str, str]:
        """Run subprocess with controlled output. All tool output is structured and timestamped."""
        if show_command or self.verbosity >= self.DEBUG:
            self.debug(f"Command: {' '.join(cmd)}")
        
        if capture_output:
            if self.verbosity >= self.DEBUG:
                # DEBUG mode: capture output, print through structured block (no raw sys.stderr writes)
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    bufsize=1,
                    encoding='utf-8',
                    errors='replace'
                )
                
                stdout_lines = []
                stderr_lines = []
                
                def read_debug_stream(stream, line_list, tag):
                    """Collect lines without writing to sys.stderr/stdout directly."""
                    for line in stream:
                        line_list.append(line)
                
                t1 = threading.Thread(target=read_debug_stream, args=(process.stdout, stdout_lines, "STDOUT"), daemon=True)
                t2 = threading.Thread(target=read_debug_stream, args=(process.stderr, stderr_lines, "STDERR"), daemon=True)
                t1.start()
                t2.start()
                
                return_code = process.wait()
                t1.join(timeout=2)
                t2.join(timeout=2)
                
                stdout_str = "".join(stdout_lines)
                stderr_str = "".join(stderr_lines)
                
                # Print structured block for STDOUT, then STDERR
                if stdout_str.strip():
                    self._debug_block_header(operation_name, "STDOUT")
                    for line in stdout_str.splitlines():
                        self._debug_tool_line(line, f"STDOUT:{operation_name}")
                    self._debug_block_footer(operation_name, "STDOUT")
                else:
                    # Still log even if nothing to print
                    if stdout_str:
                        self._log(f"[STDOUT:{operation_name}] --- begin ---")
                        for line in stdout_str.splitlines():
                            self._log(f"[STDOUT:{operation_name}] {line}")
                        self._log(f"[STDOUT:{operation_name}] --- end ---")

                if stderr_str.strip():
                    self._debug_block_header(operation_name, "STDERR")
                    for line in stderr_str.splitlines():
                        self._debug_tool_line(line, f"STDERR:{operation_name}")
                    self._debug_block_footer(operation_name, "STDERR")
                else:
                    if stderr_str:
                        self._log(f"[STDERR:{operation_name}] --- begin ---")
                        for line in stderr_str.splitlines():
                            self._log(f"[STDERR:{operation_name}] {line}")
                        self._log(f"[STDERR:{operation_name}] --- end ---")
                
                return return_code, stdout_str, stderr_str
            else:
                # Non-debug: fully suppress terminal output, log everything
                result = subprocess.run(cmd, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE, text=True, encoding='utf-8')
                if result.stdout:
                    self._log(f"[STDOUT:{operation_name}] --- begin ---")
                    for line in result.stdout.splitlines():
                        self._log(f"[STDOUT:{operation_name}] {line}")
                    self._log(f"[STDOUT:{operation_name}] --- end ---")
                if result.stderr:
                    self._log(f"[STDERR:{operation_name}] --- begin ---")
                    for line in result.stderr.splitlines():
                        self._log(f"[STDERR:{operation_name}] {line}")
                    self._log(f"[STDERR:{operation_name}] --- end ---")
                return result.returncode, result.stdout, result.stderr
        else:
            # capture_output=False: always capture so we can log/display cleanly
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    text=True, encoding='utf-8', errors='replace')
            stdout_str = result.stdout or ""
            stderr_str = result.stderr or ""
            if stdout_str.strip():
                self._debug_block_header(operation_name, "STDOUT")
                for line in stdout_str.splitlines():
                    self._debug_tool_line(line, f"STDOUT:{operation_name}")
                self._debug_block_footer(operation_name, "STDOUT")
            elif stdout_str:
                self._log(f"[STDOUT:{operation_name}] --- begin ---")
                for line in stdout_str.splitlines():
                    self._log(f"[STDOUT:{operation_name}] {line}")
                self._log(f"[STDOUT:{operation_name}] --- end ---")
            if stderr_str.strip():
                self._debug_block_header(operation_name, "STDERR")
                for line in stderr_str.splitlines():
                    self._debug_tool_line(line, f"STDERR:{operation_name}")
                self._debug_block_footer(operation_name, "STDERR")
            elif stderr_str:
                self._log(f"[STDERR:{operation_name}] --- begin ---")
                for line in stderr_str.splitlines():
                    self._log(f"[STDERR:{operation_name}] {line}")
                self._log(f"[STDERR:{operation_name}] --- end ---")
            return result.returncode, stdout_str, stderr_str
            
    def run_subprocess_with_progress(
        self,
        cmd: List[str],
        operation_name: str,
        total: int,
        unit: str = "items",
        parser_func: Optional[Callable[[str], Optional[int]]] = None,
        show_command: bool = False
    ) -> int:
        """
        Run subprocess with real-time progress tracking.
        In debug mode: shows structured tool output block instead of progress bar.
        """
        if show_command or self.verbosity >= self.DEBUG:
            self.debug(f"Command: {' '.join(cmd)}")
        
        # In debug mode, show structured output instead of a progress bar
        if self.verbosity >= self.DEBUG:
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
                encoding='utf-8',
                errors='replace'
            )
            
            stdout_lines = []
            stderr_lines = []

            def collect_stream(stream, line_list):
                for line in stream:
                    line_list.append(line)

            t1 = threading.Thread(target=collect_stream, args=(process.stdout, stdout_lines), daemon=True)
            t2 = threading.Thread(target=collect_stream, args=(process.stderr, stderr_lines), daemon=True)
            t1.start()
            t2.start()
            return_code = process.wait()
            t1.join(timeout=2)
            t2.join(timeout=2)

            stdout_str = "".join(stdout_lines)
            stderr_str = "".join(stderr_lines)

            if stdout_str.strip():
                self._debug_block_header(operation_name, "STDOUT")
                for line in stdout_str.splitlines():
                    self._debug_tool_line(line, f"STDOUT:{operation_name}")
                self._debug_block_footer(operation_name, "STDOUT")
            else:
                if stdout_str:
                    self._log(f"[STDOUT:{operation_name}] --- begin ---")
                    for line in stdout_str.splitlines():
                        self._log(f"[STDOUT:{operation_name}] {line}")
                    self._log(f"[STDOUT:{operation_name}] --- end ---")

            if stderr_str.strip():
                self._debug_block_header(operation_name, "STDERR")
                for line in stderr_str.splitlines():
                    self._debug_tool_line(line, f"STDERR:{operation_name}")
                self._debug_block_footer(operation_name, "STDERR")
            else:
                if stderr_str:
                    self._log(f"[STDERR:{operation_name}] --- begin ---")
                    for line in stderr_str.splitlines():
                        self._log(f"[STDERR:{operation_name}] {line}")
                    self._log(f"[STDERR:{operation_name}] --- end ---")

            return return_code

        # Standard/verbose mode: use progress bar as before
        progress = None
        if self.verbosity >= self.STANDARD:
            progress = ProgressBar(total, operation_name, unit)
        
        # Start subprocess with pipes
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            encoding='utf-8',
            errors='replace'
        )
        
        # Read output in real-time and update progress
        stdout_lines = []
        stderr_lines = []
        
        def read_stream(stream, line_list, is_stderr=False):
            """Read stream line by line and update progress"""
            tag = "STDERR" if is_stderr else "STDOUT"
            try:
                for line in stream:
                    line_list.append(line)
                    # Log with timestamp only — no terminal print in standard mode
                    self._log(f"[{tag}:{operation_name}] {line.rstrip()}")
                    # Parse progress if parser provided and progress bar exists
                    if parser_func and progress:
                        try:
                            parsed_value = parser_func(line)
                            if parsed_value is not None:
                                progress.update_to(parsed_value)
                        except Exception as e:
                            self.debug(f"Progress parsing error: {e}")
            except Exception as e:
                self.debug(f"Stream reading error: {e}")
        
        stdout_thread = threading.Thread(
            target=read_stream, args=(process.stdout, stdout_lines, False), daemon=True
        )
        stderr_thread = threading.Thread(
            target=read_stream, args=(process.stderr, stderr_lines, True), daemon=True
        )
        
        stdout_thread.start()
        stderr_thread.start()
        
        return_code = process.wait()
        stdout_thread.join(timeout=2)
        stderr_thread.join(timeout=2)
        
        if progress and return_code == 0:
            progress.update_to(total)
        
        return return_code

    # ========================================================================
    # UTILITY FUNCTIONS
    # ========================================================================
    
    def _format_time(self, seconds: float) -> str:
        """Format time duration"""
        if seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            minutes = int(seconds // 60)
            secs = seconds % 60
            return f"{minutes}m {secs:.0f}s"
        else:
            hours = int(seconds // 3600)
            minutes = int((seconds % 3600) // 60)
            return f"{hours}h {minutes}m"

    def stage_complete(self, title: Optional[str] = None):
        """Display stage completion message"""
        stage = title or self.current_stage
        if stage and stage in self.stage_times:
            elapsed = time.time() - self.stage_times[stage]
            if self.verbosity >= self.STANDARD:
                self._print(f"\n     {Colors.GREEN}✓ {stage} complete{Colors.END} ({self._format_time(elapsed)})")

    def file_list(self, title: str, files: List[str]):
        """Display file list"""
        if self.verbosity >= self.VERBOSE:
            self._print(f"\n{title}:")
            for file in files:
                self._print(f"     • {file}")

    def debug(self, message: str):
        """Print debug message"""
        if self.verbosity >= self.DEBUG:
            self._print(f"     {Colors.DIM}[DEBUG] {message}{Colors.END}")

    def _debug_tool_line(self, line: str, tag: str):
        """
        Print a single tool output line in debug mode.
        Goes to terminal as a structured dim line — never writes raw to sys.stderr/stdout.
        Also logs with timestamp via _log().
        """
        stripped = line.rstrip()
        # Print to terminal only — clean, indented, dim
        if self.verbosity >= self.DEBUG:
            print(f"  {Colors.DIM}│ {stripped}{Colors.END}")
        # Always log with timestamp (separate from terminal print)
        self._log(f"[{tag}] {stripped}")

    def _debug_block_header(self, operation_name: str, tag: str):
        """Print a clean visual block header for tool output in debug mode."""
        if self.verbosity >= self.DEBUG:
            # Stop + clear the spinner line before printing (prevents interleaving)
            if self._active_spinner and self._active_spinner.spinning:
                self._active_spinner.spinning = False
                if self._active_spinner.thread:
                    self._active_spinner.thread.join(timeout=0.3)
                # Clear the spinner line
                sys.stdout.write('\r\033[K')
                sys.stdout.flush()
            label = f" {tag}: {operation_name} "
            bar = "─" * max(0, (73 - len(label)) // 2)
            print(f"\n  {Colors.DIM}{bar}{label}{bar}{Colors.END}")
        self._log(f"[{tag}:{operation_name}] --- begin ---")

    def _debug_block_footer(self, operation_name: str, tag: str):
        """Print a clean visual block footer for tool output in debug mode."""
        if self.verbosity >= self.DEBUG:
            print(f"  {Colors.DIM}{'─' * 73}{Colors.END}\n")
        self._log(f"[{tag}:{operation_name}] --- end ---")

    def result(self, metrics: Dict[str, Any], indent: int = 1):
        """Display result metrics"""
        if self.verbosity >= self.STANDARD:
            prefix = "        " if indent == 2 else "     "
            for key, value in metrics.items():
                self._print(f"{prefix}→ {key}: {value}")


# ========================================================================
# GLOBAL FORMATTER INSTANCE
# ========================================================================
_formatter_instance = None

def get_formatter() -> OutputFormatter:
    global _formatter_instance
    if _formatter_instance is None:
        _formatter_instance = OutputFormatter()
    return _formatter_instance

def set_formatter(formatter: OutputFormatter):
    global _formatter_instance
    _formatter_instance = formatter