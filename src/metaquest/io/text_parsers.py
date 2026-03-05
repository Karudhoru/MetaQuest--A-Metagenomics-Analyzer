"""
MetaQuest Text Parsing Utilities
=================================
Functions for extracting information from text descriptions and headers.
"""

import re
from typing import Optional


def extract_organism_name(description: str) -> Optional[str]:
    """
    Extract organism name from protein description.
    
    Supports multiple annotation formats:
    - SwissProt: OS=Escherichia coli OX=562
    - NCBI: beta-lactamase [Escherichia coli]
    - GenBank: organism="Escherichia coli"
    - RefSeq: WP_123456.1 [Escherichia coli]
    
    Filters out database tags, accession numbers, and metadata.
    
    Args:
        description: Protein description string from BLAST/DIAMOND
        
    Returns:
        Normalized organism name (Genus species) or None if not found
        
    Examples:
        >>> extract_organism_name("Beta-lactamase OS=Escherichia coli OX=562")
        'Escherichia coli'
        
        >>> extract_organism_name("toxin [Staphylococcus aureus str. MRSA252]")
        'Staphylococcus aureus'
        
        >>> extract_organism_name("hypothetical protein [VFclass:VF0001]")
        None  # Filters out database tags
    """
    if not description or not isinstance(description, str):
        return None
    
    # Strategy 1: SwissProt OS= format (highest priority)
    os_match = re.search(
        r'OS=([^=]+?)(?:\s+(?:OX|GN|PE|SV)=|\s*$)', 
        description, 
        re.IGNORECASE
    )
    if os_match:
        return _normalize_species_name(os_match.group(1).strip())
    
    # Strategy 2: NCBI [Organism] format - with filtering
    bracket_matches = re.findall(r'\[([^\]]+)\]', description)
    if bracket_matches:
        # Filter out non-organism brackets
        skip_patterns = [
            'vfclass:', 'card:', 'aro:', 'refseq:', 'genbank:', 'gb|', 'ref|',
            'ec:', 'go:', 'ko:', 'taxon:', 'id:', 'accession:', 'locus:'
        ]
        
        for match in reversed(bracket_matches):  # Check from end (organism usually last)
            match_lower = match.lower()
            
            # Skip database tags
            if any(pattern in match_lower for pattern in skip_patterns):
                continue
            
            # Skip pure accession numbers (e.g., WP_123456.1)
            if re.match(r'^[A-Z]{2,}_[\d.]+$', match):
                continue
            
            # Skip short tags
            if len(match.split()) == 1 and len(match) < 10:
                continue
            
            # Valid organism: 2+ words or single word >10 chars
            if len(match.split()) >= 2 or len(match) >= 10:
                if re.search(r'[A-Za-z]{3,}', match):
                    return _normalize_species_name(match)
    
    # Strategy 3: GenBank organism="" format
    genbank_match = re.search(r'organism="([^"]+)"', description, re.IGNORECASE)
    if genbank_match:
        return _normalize_species_name(genbank_match.group(1).strip())
    
    # Strategy 4: Trailing organism name (after whitespace)
    trailing_match = re.search(r'\s{3,}([A-Z][a-z]+\s+[a-z]+)\s*$', description)
    if trailing_match:
        return _normalize_species_name(trailing_match.group(1).strip())
    
    return None


def _normalize_species_name(organism: str) -> str:
    """
    Normalize organism to Genus species level, removing strain info.
    
    Removes:
    - Strain designations: str. K-12, ATCC 25922
    - Subspecies: subsp. enterica
    - Serovars: serovar Typhimurium
    - Serotypes: O157:H7
    - Parenthetical notes: (type strain)
    
    Args:
        organism: Raw organism string
        
    Returns:
        Normalized "Genus species" string
        
    Examples:
        >>> _normalize_species_name("Escherichia coli str. K-12")
        'Escherichia coli'
        
        >>> _normalize_species_name("Salmonella enterica subsp. enterica serovar Typhimurium")
        'Salmonella enterica'
        
        >>> _normalize_species_name("Staphylococcus aureus ATCC 25923")
        'Staphylococcus aureus'
    """
    # Remove strain-specific suffixes
    patterns_to_remove = [
        r'\s+serovar\s+\S+',     # serovar Typhimurium
        r'\s+str\.\s+\S+',       # str. K-12
        r'\s+substr\.\s+\S+',    # substr. thermophilus
        r'\s+subsp\.\s+\S+',     # subsp. enterica
        r'\s+strain\s+\S+',      # strain NCTC 8325
        r'\s+[A-Z][0-9]+:[HO][0-9]+',  # O157:H7
        r'\s+ATCC\s+\S+',        # ATCC 25922
        r'\s+DSM\s+\S+',         # DSM 12345
        r'\s+NCTC\s+\S+',        # NCTC 8325
        r'\s+\([^)]+\)',         # (anything in parentheses)
    ]
    
    for pattern in patterns_to_remove:
        organism = re.sub(pattern, '', organism, flags=re.IGNORECASE)
    
    # Clean up whitespace
    organism = ' '.join(organism.split())
    
    # Keep only Genus species (first two words)
    parts = organism.split()
    if len(parts) >= 2:
        # Handle "Candidatus Organism species"
        if parts[0].lower() == 'candidatus' and len(parts) >= 3:
            return f"{parts[0]} {parts[1]} {parts[2]}"
        return f"{parts[0]} {parts[1]}"
    elif len(parts) == 1:
        return parts[0]  # Genus only
    
    return organism.strip()
