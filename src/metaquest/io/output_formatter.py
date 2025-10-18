"""
MetaQuest Professional Output Formatter v4.0.0
==============================================

A sophisticated, multi-level output formatting system for CLI operations.
Provides clean, structured, and professional terminal output with configurable verbosity.

Features:
- Two-tier verbosity system (standard, debug)
- Smart subprocess output capture and filtering
- Progress tracking with spinners and progress bars
- Hierarchical status indicators
- Color-aware terminal output
- Time tracking and estimation
- Actionable error messages

Author: MetaQuest Development Team
"""

import sys
import time
import subprocess
import threading
import re
from pathlib import Path
from datetime import datetime, timedelta
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
    """Animated spinner for long-running operations"""
    def __init__(self, message: str = "Processing"):
        self.message = message
        self.spinning = False
        self.thread = None
        self.frames = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏']
        self.frame_idx = 0
    def _spin(self):
        """Spinner animation loop"""
        while self.spinning:
            frame = self.frames[self.frame_idx % len(self.frames)]
            sys.stdout.write(f'\r [{frame}] {self.message}...')
            sys.stdout.flush()
            self.frame_idx += 1
            time.sleep(0.1)
    def start(self):
        """Start the spinner"""
        self.spinning = True
        self.thread = threading.Thread(target=self._spin)
        self.thread.daemon = True
        self.thread.start()
    def stop(self, final_message: Optional[str] = None):
        """Stop the spinner"""
        self.spinning = False
        if self.thread:
            self.thread.join()
        if final_message:
            sys.stdout.write(f'\r {final_message}\n')
        else:
            sys.stdout.write('\r' + ' ' * (len(self.message) + 10) + '\r')
        sys.stdout.flush()

class ProgressBar:
    """Progress bar for tracking operations"""
    def __init__(self, total: int, desc: str = "", unit: str = "it", width: int = 40):
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
        """Display the progress bar"""
        percent = self.current / self.total
        filled = int(self.width * percent)
        bar = '█' * filled + '░' * (self.width - filled)

        elapsed = time.time() - self.start_time
        rate = self.current / elapsed if elapsed > 0 else 0
        remaining_s = (self.total - self.current) / rate if rate > 0 else 0
        remaining_str = self._format_time(remaining_s) if self.current < self.total else "done"

        status = (f"\r    {self.desc} {bar} {percent:3.0%}"
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

class OutputFormatter:
    """
    Professional output formatter with two-level verbosity support

    Verbosity Levels:
    - standard: Clean, essential information + operational details
    - debug: Full diagnostic output
    """

    STANDARD = 1
    DEBUG = 2

    def __init__(self, verbosity: str = 'standard', log_file: Optional[Path] = None):
        """
        Initialize the output formatter

        Args:
            verbosity: One of 'standard', 'debug'
            log_file: Optional path to write full logs
        """
        self.verbosity_map = {
            'standard': self.STANDARD,
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
        if not self.colors_enabled:
            Colors.disable()
        # Open log file if specified
        if self.log_file:
            self.log_file.parent.mkdir(parents=True, exist_ok=True)
            self.log_handle = open(self.log_file, 'w', encoding='utf-8')
        else:
            self.log_handle = None

    def __del__(self):
        """Clean up log file handle"""
        if hasattr(self, 'log_handle') and self.log_handle:
            self.log_handle.close()

    def _log(self, message: str):
        """Write to log file if enabled"""
        if self.log_handle:
            timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            self.log_handle.write(f"[{timestamp}] {message}\n")
            self.log_handle.flush()

    def _print(self, message: str, min_verbosity: int = STANDARD):
        """Internal print with verbosity check"""
        if self.verbosity >= min_verbosity:
            print(message)
        self._log(message)

    # ========================================================================
    # SECTION HEADERS
    # ========================================================================
    def banner(self, title: str, version: str, tagline: str):
        """Print application banner"""
        if self.verbosity >= self.STANDARD:
            banner = f"""
{Colors.BOLD}{Colors.CYAN}{'═' * 75}
  {title} v{version}
  {tagline}
  MetaQuest Development Team
{'═' * 75}{Colors.END}
"""
            print(banner)
            self._log(f"{title} v{version} - {tagline}")

    def step_header(self, step_num: int, total_steps: int, title: str):
        """Print major pipeline step header"""
        self.current_stage = title
        self.stage_times[title] = time.time()
        if self.verbosity >= self.STANDARD:
            header = f"\n{'═' * 70}\nSTEP {step_num}/{total_steps}: {title}\n{'═' * 70}"
            self._print(header)

    def section_header(self, title: str):
        """Print subsection header"""
        if self.verbosity >= self.STANDARD:
            self._print(f"\n{Colors.BOLD}{Colors.CYAN}[{title}]{Colors.END}")
            self._print(f"{Colors.CYAN}{'─' * (len(title) + 2)}{Colors.END}")

    def header(self, title: str):
        self.section_header(title)

    # ========================================================================
    # STATUS MESSAGES
    # ========================================================================
    def operation(self, message: str, show_in_standard: bool = True):
        """
        Print operation start message
        Returns an operation ID for tracking completion
        """
        min_level = self.STANDARD if show_in_standard else self.DEBUG
        if self.verbosity >= min_level:
            self._print(f"   [•] {message}...")
        op_id = len(self._operation_stack)
        self._operation_stack.append(message)
        return op_id

    def step(self, message: str):
        return self.operation(message, show_in_standard=True)

    def substep(self, message: str):
        if self.verbosity >= self.DEBUG:
            self._print(f"      → {message}...")

    def live_status(self, message: str):
        """
        Prints a status message that can be overwritten.
        Ensures the line is cleared to the end.
        """
        if self.verbosity >= self.STANDARD:
            # Get terminal width to clear the entire line
            terminal_width = shutil.get_terminal_size((80, 20)).columns
            # \r moves cursor to the start, message is printed, and ANSI clear ensures old text is gone
            sys.stdout.write(f"\r   [ ] {message}".ljust(terminal_width))
            sys.stdout.flush()

    def clear_line(self):
        """Clears the current line after a live_status update is complete."""
        if self.verbosity >= self.STANDARD:
            terminal_width = shutil.get_terminal_size((80, 20)).columns
            sys.stdout.write("\r".ljust(terminal_width) + "\r")
            sys.stdout.flush()

    def success(self, message: str, style: Optional[str] = None):
        if self.verbosity >= self.STANDARD:
            prefix = f"   {Colors.GREEN}✓{Colors.END} "
            if style == 'bold':
                message = f"{Colors.BOLD}{message}{Colors.END}"
            self._print(f"{prefix}{message}")

    def complete(self, message: str):
        self.success(message)

    def info(self, message: str, indent: int = 1):
        if self.verbosity >= self.STANDARD:
            prefix = "   " * indent + "→ "
            self._print(f"{prefix}{message}")

    def metric(self, key: str, value: Any, indent: int = 1):
        if self.verbosity >= self.STANDARD:
            prefix = "   " * indent + "→ "
            self._print(f"{prefix}{key}: {value}")

    def warning(self, message: str):
        if self.verbosity >= self.STANDARD:
            self._print(f"   {Colors.YELLOW}⚠{Colors.END}  {message}")

    def error(self, message: str, solutions: Optional[List[str]] = None, doc_link: Optional[str] = None):
        if self.verbosity >= self.STANDARD:
            self._print(f"\n   {Colors.RED}✗{Colors.END} {Colors.RED}{message}{Colors.END}")
            if solutions:
                self._print(f"\n   {Colors.BOLD}Solutions:{Colors.END}")
                for solution in solutions:
                    self._print(f"   • {solution}")
            if doc_link:
                self._print(f"\n   {Colors.CYAN}Documentation:{Colors.END} {doc_link}")
            print()

    def debug(self, message: str):
        if self.verbosity >= self.DEBUG:
            self._print(f"   {Colors.DIM}[DEBUG] {message}{Colors.END}")

    # ========================================================================
    # METRICS & RESULTS
    # ========================================================================
    def result(self, metrics: Dict[str, Any], indent: int = 1):
        if self.verbosity >= self.STANDARD:
            prefix = "   " * indent + "→ "
            for key, value in metrics.items():
                self._print(f"{prefix}{key}: {value}")

    def summary(self, title: str, metrics: Dict[str, Any]):
        if self.verbosity >= self.STANDARD:
            self._print(f"\n   {'─' * 66}")
            self._print(f"   {Colors.BOLD}✓ {title}{Colors.END}")
            self._print(f"   {'─' * 66}")
            for key, value in metrics.items():
                self._print(f"      • {key}: {value}")
            self._print(f"   {'─' * 66}")

    def stage_complete(self, title: Optional[str] = None):
        stage = title or self.current_stage
        if stage and stage in self.stage_times:
            elapsed = time.time() - self.stage_times[stage]
            if self.verbosity >= self.STANDARD:
                self._print(f"\n   {'─' * 66}")
                self._print(f"   {Colors.GREEN}✓ {stage} complete{Colors.END}")
                self._print(f"   Time: {self._format_time(elapsed)}")
                self._print(f"   {'─' * 66}")

    # ========================================================================
    # PROGRESS TRACKING
    # ========================================================================

    @contextmanager
    def spinner(self, message: str):
        if self.verbosity >= self.STANDARD and self.colors_enabled:
            # Stop any live status to prevent overlap
            self.clear_line() 
            spinner = Spinner(message)
            self.is_spinning = True
            spinner.start()
            try:
                yield
            finally:
                spinner.stop()
                self.is_spinning = False
        else:
            self.operation(message)
            yield

    def progress_bar(self, total: int, desc: str, unit: str = "it") -> ProgressBar:
        if self.verbosity >= self.STANDARD:
            return ProgressBar(total, desc, unit)
        else:
            class NoOpProgressBar: # A silent progress bar for non-standard modes
                def update(self, n=1): pass
                def update_to(self, v): pass
                def __enter__(self): return self
                def __exit__(self, *args): pass
            return NoOpProgressBar()

    # ========================================================================
    # SUBPROCESS EXECUTION
    # ========================================================================
    def run_subprocess(self, cmd: List[str], operation_name: str, capture_output: bool = True, show_command: bool = False) -> Tuple[int, str, str]:
        if show_command:
            self.debug(f"Command: {' '.join(cmd)}")
        
        if capture_output:
            # Only start a spinner if one isn't already active and we are in standard mode.
            if self.verbosity >= self.STANDARD and not self.is_spinning:
                with self.spinner(operation_name):
                    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding='utf-8')
            else:
                # Otherwise, run silently (the outer spinner is still visible to the user).
                result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding='utf-8')
            
            if result.stdout:
                self._log(f"STDOUT for '{operation_name}':\n{result.stdout}")
            if result.stderr:
                self._log(f"STDERR for '{operation_name}':\n{result.stderr}")
            return result.returncode, result.stdout, result.stderr
        else:
            # Logic for non-captured output remains the same
            if self.verbosity >= self.DEBUG:
                result = subprocess.run(cmd)
                return result.returncode, "", ""
            else:
                result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if result.stdout:
                    self._log(f"STDOUT:\n{result.stdout.decode(errors='ignore')}")
                if result.stderr:
                    self._log(f"STDERR:\n{result.stderr.decode(errors='ignore')}")
                return result.returncode, "", ""
            
    # Replace the existing run_subprocess_with_progress with this new version

    def run_subprocess_with_progress(
    self,
    cmd: List[str],
    operation_name: str,
    total: int,
    unit: str,
    parser_func: Callable[[str], Optional[int]]
) -> int:
        """
        Runs a subprocess and drives a progress bar by parsing its live output.

        Args:
            cmd: Command to run.
            operation_name: High-level description of the task.
            total: The total number of items for the progress bar (e.g., number of proteins).
            unit: The unit for the progress bar (e.g., "proteins").
            parser_func: A function that takes a line of output and returns an integer
                        representing the *absolute* progress (e.g., 20000).

        Returns:
            The final return code of the subprocess.
        """
        self.operation(operation_name)
        self.debug(f"Command: {' '.join(cmd)}")

        # Use unbuffered mode and line buffering
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding='utf-8',
            bufsize=1,  # Line buffered
            universal_newlines=True
        )

        last_progress = 0
        update_count = 0
        
        with self.progress_bar(total=total, desc=f"      {operation_name}", unit=unit) as pbar:
            try:
                for line in iter(process.stdout.readline, ''):
                    if not line:  # Empty line means EOF
                        break
                        
                    # Always log the line
                    self._log(line.rstrip())
                    
                    # In debug mode, print the tool's raw output
                    if self.verbosity >= self.DEBUG:
                        print(f"    {Colors.DIM}{line.rstrip()}{Colors.END}")

                    # Try to parse progress
                    progress_value = parser_func(line)
                    if progress_value is not None and progress_value > last_progress:
                        # Update the progress bar to the absolute value parsed
                        pbar.update_to(progress_value)
                        last_progress = progress_value
                        update_count += 1
            
            except KeyboardInterrupt:
                self.warning("Interrupted by user, terminating process...")
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                raise
            
            finally:
                # Ensure process is cleaned up
                if process.poll() is None:
                    process.terminate()
        
        return_code = process.wait()
        
        # Debug info about progress updates
        if update_count == 0:
            self.debug(f"Warning: No progress updates detected during {operation_name}")
            self.debug("Consider checking if the tool outputs progress information")
        else:
            self.debug(f"Progress bar updated {update_count} times")
        
        if return_code != 0:
            self.warning(f"'{operation_name}' finished with a non-zero exit code: {return_code}")
        
        return return_code

    # ========================================================================
    # UTILITY FUNCTIONS
    # ========================================================================
    def _format_time(self, seconds: float) -> str:
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

    def get_elapsed_time(self) -> str:
        elapsed = time.time() - self.start_time
        return self._format_time(elapsed)

    def file_list(self, title: str, files: List[str]):
        if self.verbosity >= self.STANDARD:
            self._print(f"\n{title}:")
            for file in files:
                self._print(f"   • {file}")

# ========================================================================
# OUTPUT PARSERS FOR COMMON TOOLS
# ========================================================================
class ToolOutputParsers:
    """Collection of output parsers for common bioinformatics tools"""
    @staticmethod
    def parse_kraken2(line: str) -> Optional[Dict[str, Any]]:
        if "sequences" in line and "processed" in line:
            match = re.search(r'(\d+)\s+sequences.*in\s+([\d.]+)s', line)
            if match:
                return {'sequences_processed': int(match.group(1)), 'time': float(match.group(2))}
        if "classified" in line:
            match = re.search(r'(\d+)\s+sequences classified.*\(([\d.]+)%\)', line)
            if match:
                return {'classified': int(match.group(1)), 'classified_pct': float(match.group(2))}
        return None
    @staticmethod
    def parse_bracken(line: str) -> Optional[Dict[str, Any]]:
        if "Number of species in sample:" in line:
            match = re.search(r'Number of species in sample:\s*(\d+)', line)
            if match:
                return {'species_detected': int(match.group(1))}
        if "Total reads in sample:" in line:
            match = re.search(r'Total reads in sample:\s*(\d+)', line)
            if match:
                return {'total_reads': int(match.group(1))}
        return None
    @staticmethod
    def parse_diamond(line: str) -> Optional[Dict[str, Any]]:
        if "Reported" in line and "alignments" in line:
            match = re.search(r'Reported\s+(\d+)\s+pairwise alignments', line)
            if match:
                return {'alignments': int(match.group(1))}
        if "queries aligned" in line:
            match = re.search(r'(\d+)\s+queries aligned', line)
            if match:
                return {'queries_aligned': int(match.group(1))}
        return None
    @staticmethod
    def parse_prokka(line: str) -> Optional[Dict[str, Any]]:
        if "annotation:" in line.lower():
            pass
        return None

# ========================================================================
# GLOBAL FORMATTER INSTANCE (Singleton pattern)
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
