#!/usr/bin/env python3

"""
MetaQuest Visualization Module 
Base Visualizer Class
"""

import json
from pathlib import Path

class BaseVisualizer:
    """Base class for all visualizers with common functionality"""
    
    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def save_plot(self, fig, filename: str) -> str:
        """Save plotly figure to HTML file"""
        filepath = self.output_dir / filename
        fig.write_html(filepath)
        return str(filepath)
    
    def load_json_report(self, report_path: Path) -> dict:
        """Load and validate JSON report file"""
        try:
            if not report_path.exists():
                print(f"⚠️ Report file not found: {report_path}")
                return {}
            
            if report_path.stat().st_size == 0:
                print(f"⚠️ Report file is empty: {report_path}")
                return {}
            
            with open(report_path, 'r') as f:
                data = json.load(f)
            
            return data
        except json.JSONDecodeError as e:
            print(f"⚠️ Invalid JSON in report file: {e}")
            return {}
        except Exception as e:
            print(f"⚠️ Error loading report: {e}")
            return {}
