#!/usr/bin/env python3

"""
MetaQuest Visualization Module 
Base Visualizer Class
*** REFACTORED to use global OutputFormatter ***
"""

import json
from pathlib import Path
from ..io.output_formatter import get_formatter # <-- IMPORT FORMATTER

class BaseVisualizer:
    """Base class for all visualizers with common functionality"""
    
    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.formatter = get_formatter() # <-- ADD FORMATTER
    
    def save_plot(self, fig, filename: str) -> str:
        """Save plotly figure to HTML file"""
        filepath = self.output_dir / filename
        fig.write_html(filepath)
        # Note: The *calling* function should log success, e.g.:
        # self.formatter.success(f"Plot saved: {filename}")
        return str(filepath)
    
    def load_json_report(self, report_path: Path) -> dict:
        """Load and validate JSON report file"""
        try:
            if not report_path.exists():
                self.formatter.warning(f"Report file not found: {report_path}") # <-- CHANGE
                return {}
            
            if report_path.stat().st_size == 0:
                self.formatter.warning(f"Report file is empty: {report_path}") # <-- CHANGE
                return {}
            
            with open(report_path, 'r') as f:
                data = json.load(f)
            
            return data
        except json.JSONDecodeError as e:
            self.formatter.error(f"Invalid JSON in report file: {e}") # <-- CHANGE
            return {}
        except Exception as e:
            self.formatter.error(f"Error loading report: {e}") # <-- CHANGE
            return {}