"""
MetaQuest Base Reporter - Professional Edition v5.0.0
=======================================================
Foundation class for all MetaQuest report generation.

ENHANCEMENTS v5.0.0:
  ✓ Unified version string across all reporters
  ✓ Enhanced professional formatting helpers
  ✓ Publication-ready table generators
  ✓ Consistent decimal precision handling
  ✓ JSON/CSV export utilities
  ✓ File I/O helpers
  
Author: MetaQuest Development Team
License: MIT
"""

from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, List, Any, Union
import json
import csv


class BaseReporter:
    """
    Base class for all MetaQuest reporters.
    Provides common formatting, file I/O, and data export utilities.
    """
    
    VERSION = "5.0.0"  # Unified version across all reporters
    
    def __init__(self, output_dir: Path):
        """
        Initialize base reporter.
        
        Args:
            output_dir: Directory for output files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    # ========================================================================
    # FILE I/O UTILITIES
    # ========================================================================
    
    def save_report(self, content: str, filename: str) -> Path:
        """
        Save text report to file.
        
        Args:
            content: Report content string
            filename: Output filename
        
        Returns:
            Path to saved file
        """
        filepath = self.output_dir / filename
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
        return filepath
    
    def save_json(self, data: Union[Dict, List], filename: str, indent: int = 2) -> Path:
        """
        Save data as JSON file with pretty formatting.
        
        Args:
            data: Dictionary or list to save
            filename: Output filename (will add .json if missing)
            indent: JSON indentation level
        
        Returns:
            Path to saved file
        """
        if not str(filename).endswith('.json'):
            filename += '.json'
        
        filepath = self.output_dir / filename
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=indent, ensure_ascii=False)
        return filepath
    
    def load_json(self, filename: str) -> Union[Dict, List]:
        """
        Load JSON file from output directory.
        
        Args:
            filename: JSON filename
        
        Returns:
            Loaded data structure
        """
        filepath = self.output_dir / filename
        if not filepath.exists():
            raise FileNotFoundError(f"JSON file not found: {filepath}")
        
        with open(filepath, 'r', encoding='utf-8') as f:
            return json.load(f)
    
    def save_csv(self, data: List[Dict], filename: str, headers: Optional[List[str]] = None) -> Path:
        """
        Save data as CSV file.
        
        Args:
            data: List of dictionaries (one per row)
            filename: Output filename (will add .csv if missing)
            headers: Column headers (auto-detected from first row if None)
        
        Returns:
            Path to saved file
        """
        if not str(filename).endswith('.csv'):
            filename += '.csv'
        
        if not data:
            raise ValueError("Cannot save empty data to CSV")
        
        filepath = self.output_dir / filename
        
        # Auto-detect headers from first row if not provided
        if headers is None:
            headers = list(data[0].keys())
        
        with open(filepath, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=headers)
            writer.writeheader()
            writer.writerows(data)
        
        return filepath
    
    def load_csv(self, filename: str) -> List[Dict]:
        """
        Load CSV file as list of dictionaries.
        
        Args:
            filename: CSV filename
        
        Returns:
            List of row dictionaries
        """
        filepath = self.output_dir / filename
        if not filepath.exists():
            raise FileNotFoundError(f"CSV file not found: {filepath}")
        
        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            return list(reader)
    
    # ========================================================================
    # FORMATTING UTILITIES
    # ========================================================================
    
    @staticmethod
    def format_header(title: str, width: int = 78, style: str = "double") -> str:
        """
        Generate professional report header with consistent styling.
        
        Args:
            title: Header title text
            width: Total width of header (default: 78 for compatibility)
            style: Border style ('double', 'single', 'section', 'minimal')
        
        Returns:
            Formatted header string
        """
        if style == "double":
            top_border = "═" * width
            bottom_border = "═" * width
        elif style == "single":
            top_border = "─" * width
            bottom_border = "─" * width
        elif style == "section":
            top_border = ""
            bottom_border = "─" * width
        elif style == "minimal":
            return f"\n{title}\n{'─' * len(title)}"
        else:
            top_border = bottom_border = "=" * width
        
        # Center title with proper padding
        padding = (width - len(title)) // 2
        centered_title = " " * padding + title + " " * (width - len(title) - padding)
        
        if style == "section":
            return f"{centered_title}\n{bottom_border}"
        else:
            return f"{top_border}\n{centered_title}\n{bottom_border}"
    
    @staticmethod
    def format_section(title: str, level: int = 1, width: int = 78) -> str:
        """
        Generate hierarchical section headers.
        
        Args:
            title: Section title
            level: Section level (1=main, 2=subsection, 3=sub-subsection)
            width: Total width
        
        Returns:
            Formatted section header
        """
        if level == 1:
            return f"\n{'═' * width}\n{title.upper()}\n{'═' * width}\n"
        elif level == 2:
            return f"\n{'─' * width}\n{title}\n{'─' * width}\n"
        else:
            return f"\n{title}\n{'-' * len(title)}\n"
    
    @staticmethod
    def format_metadata(sample_id: str, analysis_date: Optional[str] = None, **kwargs) -> str:
        """
        Generate standardized metadata section for reports.
        
        Args:
            sample_id: Sample identifier
            analysis_date: ISO format datetime (auto-generated if None)
            **kwargs: Additional metadata key-value pairs
        
        Returns:
            Formatted metadata string
        """
        if analysis_date is None:
            analysis_date = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
        
        lines = [
            f"Sample ID:          {sample_id}",
            f"Analysis Date:      {analysis_date}",
        ]
        
        # Add additional metadata in order
        for key, value in kwargs.items():
            # Format key: convert snake_case to Title Case
            formatted_key = key.replace('_', ' ').title()
            # Right-pad key to 20 chars for alignment
            lines.append(f"{formatted_key:<20}{value}")
        
        lines.append(f"Pipeline Version:   MetaQuest v{BaseReporter.VERSION}")
        
        return "\n".join(lines)
    
    @staticmethod
    def format_table(headers: list, rows: list, alignments: Optional[list] = None, 
                     col_widths: Optional[list] = None, show_divider: bool = True) -> str:
        """
        Generate professional publication-ready tables with precise alignment.
        
        Args:
            headers: List of column headers
            rows: List of row data (each row is a list matching headers)
            alignments: List of alignment specs ('left', 'right', 'center')
            col_widths: List of column widths (auto-calculated if None)
            show_divider: Add divider line after header
        
        Returns:
            Formatted table string
        
        Example:
            >>> headers = ['Species', 'Abundance', 'Risk']
            >>> rows = [['E. coli', '75.48%', 'High']]
            >>> print(format_table(headers, rows, alignments=['left', 'right', 'left']))
        """
        if not rows:
            return "(No data available)"
        
        num_cols = len(headers)
        
        # Default alignments: left for text, right for numbers
        if alignments is None:
            alignments = ['left'] * num_cols
        
        # Calculate column widths
        if col_widths is None:
            col_widths = [len(str(h)) for h in headers]
            for row in rows:
                for i, cell in enumerate(row):
                    col_widths[i] = max(col_widths[i], len(str(cell)))
        
        # Build format strings for each column
        formats = []
        for i, align in enumerate(alignments):
            width = col_widths[i]
            if align == 'right':
                formats.append(f"{{:>{width}}}")
            elif align == 'center':
                formats.append(f"{{:^{width}}}")
            else:  # left
                formats.append(f"{{:<{width}}}")
        
        # Generate header
        header_line = "  ".join(formats[i].format(str(headers[i])) for i in range(num_cols))
        
        lines = [header_line]
        
        # Add divider
        if show_divider:
            total_width = sum(col_widths) + 2 * (num_cols - 1)  # 2 spaces between cols
            lines.append("─" * total_width)
        
        # Add data rows
        for row in rows:
            row_line = "  ".join(formats[i].format(str(row[i])) for i in range(num_cols))
            lines.append(row_line)
        
        return "\n".join(lines)
    
    @staticmethod
    def format_key_value(pairs: List[tuple], key_width: int = 30,
                        value_align: str = 'left', indent: int = 0) -> str:
        """
        Format key-value pairs in aligned columns.

        Args:
            pairs: List of (key, value) tuples
            key_width: Width of key column
            value_align: Alignment for values ('left' or 'right')
            indent: Number of spaces to prepend to each line (default: 0)

        Returns:
            Formatted string
        """
        prefix = " " * indent
        lines = []
        for key, value in pairs:
            if value_align == 'right':
                lines.append(f"{prefix}{key:<{key_width}}{str(value):>20}")
            else:
                lines.append(f"{prefix}{key:<{key_width}}{value}")
        return "\n".join(lines)
    
    # ========================================================================
    # NUMBER FORMATTING
    # ========================================================================
    
    @staticmethod
    def format_number(value: float, decimals: int = 2, percentage: bool = False, 
                     suffix: str = "") -> str:
        """
        Format numbers with consistent precision for publication.
        
        Args:
            value: Numeric value to format
            decimals: Number of decimal places
            percentage: Add '%' suffix if True
            suffix: Additional suffix (e.g., ' bp', ' reads')
        
        Returns:
            Formatted number string
        """
        formatted = f"{value:.{decimals}f}"
        if percentage:
            formatted += "%"
        if suffix:
            formatted += suffix
        return formatted
    
    @staticmethod
    def format_large_number(value: int, use_commas: bool = True) -> str:
        """
        Format large integers with thousand separators.
        
        Args:
            value: Integer value
            use_commas: Add comma separators
        
        Returns:
            Formatted number string (e.g., '1,234,567')
        """
        if use_commas:
            return f"{value:,}"
        return str(value)
    
    @staticmethod
    def format_scientific(value: float, precision: int = 2) -> str:
        """
        Format number in scientific notation.
        
        Args:
            value: Numeric value
            precision: Decimal places
        
        Returns:
            Scientific notation string (e.g., '1.23e-05')
        """
        return f"{value:.{precision}e}"
    
    # ========================================================================
    # VISUAL ELEMENTS
    # ========================================================================
    
    @staticmethod
    def create_bar(value: float, max_value: float, width: int = 50, 
                   char: str = "█", show_value: bool = True) -> str:
        """
        Create ASCII progress/abundance bars for visualizing proportions.
        
        Args:
            value: Current value
            max_value: Maximum value (100% equivalent)
            width: Total bar width in characters
            char: Character to use for filled portion
            show_value: Append numeric value to bar
        
        Returns:
            Formatted bar string
        """
        if max_value == 0:
            filled = 0
        else:
            proportion = min(value / max_value, 1.0)
            filled = int(proportion * width)
        
        bar = char * filled + " " * (width - filled)
        
        if show_value:
            percentage = (value / max_value * 100) if max_value > 0 else 0
            return f"{bar} {percentage:.1f}%"
        return bar
    
    @staticmethod
    def create_bullet_list(items: List[str], bullet: str = "•", indent: int = 2) -> str:
        """
        Create formatted bullet list.
        
        Args:
            items: List of items
            bullet: Bullet character
            indent: Indentation spaces
        
        Returns:
            Formatted list string
        """
        indent_str = " " * indent
        return "\n".join(f"{indent_str}{bullet} {item}" for item in items)
    
    @staticmethod
    def create_numbered_list(items: List[str], indent: int = 2) -> str:
        """
        Create formatted numbered list.
        
        Args:
            items: List of items
            indent: Indentation spaces
        
        Returns:
            Formatted list string
        """
        indent_str = " " * indent
        return "\n".join(f"{indent_str}{i}. {item}" for i, item in enumerate(items, 1))
    
    # ========================================================================
    # REPORT STRUCTURE
    # ========================================================================
    
    @staticmethod
    def format_footer(include_timestamp: bool = True, 
                     contact_email: str = "metaquest-support@example.org") -> str:
        """
        Generate standardized report footer.
        
        Args:
            include_timestamp: Add generation timestamp
            contact_email: Support contact email
        
        Returns:
            Formatted footer string
        """
        lines = []
        
        lines.append("═" * 78)
        
        if include_timestamp:
            timestamp = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
            lines.append(f"Report generated: {timestamp}")
        
        lines.append(f"MetaQuest version: v{BaseReporter.VERSION}")
        lines.append(f"For questions: {contact_email}")
        lines.append("")
        lines.append("═" * 78)
        
        return "\n".join(lines)
    
    def format_single_kv(self, key: str, value: Any = None, indent: int = 1, key_width: int = 20) -> str:
        """
        Format a single key-value pair with indentation and alignment.

        NOTE: Renamed from format_key_value to avoid overriding the @staticmethod
        of the same name that accepts a list of (key, value) tuples.

        Args:
            key: The key/label
            value: The value to display (optional)
            indent: Indentation level (1 = "   →", 2 = "      ", etc.)
            key_width: Width for key alignment (default: 20)

        Returns:
            Formatted string
        """
        # Determine prefix based on indent level
        if indent == 1:
            prefix = "   → "
        elif indent == 2:
            prefix = "      "
        else:
            prefix = "   " * indent

        # If no value provided, just return the key with prefix
        if value is None:
            return f"{prefix}{key}"

        # Align the key to the specified width
        aligned_key = f"{key}:".ljust(key_width)

        return f"{prefix}{aligned_key} {value}"



    
    def get_version(self) -> str:
        """Get current pipeline version."""
        return self.VERSION
