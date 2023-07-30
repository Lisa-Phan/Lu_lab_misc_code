import re

def correct_line_format(line):
    # Use regular expressions to split the line into different parts
    match = re.match(r'([A-Z]+)(\d+)', line)
    
    if match:
        prefix = match.group(1)
        number = match.group(2)
        
        # Insert spaces between the prefix and number
        line = f"{prefix} {number}"
    
    # Replace any dashes (-) with spaces
    line = line.replace('-', ' ')

    # Return the corrected line
    return line

# Example line for demonstration
line = "HETATM14132 "
line2 = "HETATM14132 "

# Correct the formatting of the line
corrected_line = correct_line_format(line)

print(corrected_line)
print(correct_line_format(line2))