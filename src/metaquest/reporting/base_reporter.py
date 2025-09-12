#!/usr/bin/env python3
"""
MetaQuest Professional Reporting Module.
"""
import json
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List



class BaseReportGenerator(ABC):
    """Base class for all MetaQuest report generators with professional formatting"""
    
    def __init__(self, output_dir: Path, analysis_type: str):
        self.output_dir = Path(output_dir)
        self.analysis_type = analysis_type
        self.timestamp = datetime.now()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.version = "3.5.0"
    
    @abstractmethod
    def generate_report(self, data: Any, **kwargs) -> Dict[str, str]:
        """Generate the specific report type"""
        pass
    
    def _create_professional_header(self, title: str, pipeline_info: str = None) -> List[str]:
        """Create professional scientific report header"""
        header_line = "=" * len(title)
        
        header = [
            header_line,
            title,
            header_line,
            "",
            f"Date: {self.timestamp.strftime('%Y-%m-%d %H:%M:%S')}",
            f"Pipeline: {pipeline_info or 'MetaQuest Taxonomic Classification'}",
            f"Analysis ID: MQ_TAX_{self.timestamp.strftime('%Y%m%d_%H%M%S')}",
            f"Version: MetaQuest v{self.version}",
            ""
        ]
        
        return header
    
    def _save_text_report(self, content: List[str], filename: str) -> str:
        """Save text report with standard formatting"""
        filepath = self.output_dir / filename
        with open(filepath, 'w') as f:
            f.write('\n'.join(content))
        return str(filepath)
    
    def _save_json_report(self, data: Dict, filename: str) -> str:
        """Save JSON report with metadata"""
        report_data = {
            'metadata': {
                'analysis_type': self.analysis_type,
                'generation_timestamp': self.timestamp.isoformat(),
                'pipeline_version': f'MetaQuest v{self.version}'
            },
            'data': data
        }
        filepath = self.output_dir / filename
        with open(filepath, 'w') as f:
            json.dump(report_data, f, indent=2)
        return str(filepath)
