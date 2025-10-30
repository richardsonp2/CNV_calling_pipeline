#!/usr/bin/env python3
"""
Parse YAML config and output as bash variables
Usage: eval $(python parse_config.py penncnv_config.yaml)
"""

import sys
import yaml

def flatten_dict(d, parent_key='', sep='_'):
    """Flatten nested dictionary into bash variable names"""
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}".upper() if parent_key else k.upper()
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        elif isinstance(v, bool):
            items.append((new_key, 'true' if v else 'false'))
        else:
            items.append((new_key, str(v)))
    return dict(items)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: parse_config.py <config.yaml>", file=sys.stderr)
        sys.exit(1)
    
    with open(sys.argv[1], 'r') as f:
        config = yaml.safe_load(f)
    
    flat_config = flatten_dict(config)
    
    # Output as bash variable assignments
    for key, value in flat_config.items():
        # Escape special characters for bash
        safe_value = str(value).replace('"', '\\"')
        print(f'export CFG_{key}="{safe_value}"')
