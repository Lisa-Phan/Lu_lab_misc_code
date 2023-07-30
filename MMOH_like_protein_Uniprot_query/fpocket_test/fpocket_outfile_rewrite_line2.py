import re

def correct_HETAM_format(line):
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

def correct_number_format(line):
    # Use regular expressions to replace the minus sign (-) with a space
    line = re.sub(r'(\S)-(\S)', r'\1 -\2', line)

    # Return the corrected line
    return line
import re

def line_formatting(line):
    # Add a space between the 'HETATM' prefix and the number
    line = re.sub(r'(HETATM)(\d+)', r'\1 \2', line)

    # Replace the minus sign (-) with a space
    line = re.sub(r'(\S)-(\S)', r'\1 -\2', line)

    # Return the corrected line
    return line


# Example line for demonstration
line = "HETATM13988 APOL STP C  1 133.371-106.082 -50.464"
line2 = "HETATM13988 APOL STP C 133.371 -106.082 -50.464"
line3 = "HETATM13988 APOL STP C 133.371 -106.082-50.464"
line4 = "HETATM13988 APOL STP C 133.371-106.082-50.464"
line5 = "HETATM13988 APOL STP C   1     167.215-125.930 -35.809  0.00  0.00          Ve" 
line6 = "HETATM 13988 APOL STP C   1     167.215 -125.930 -35.809  0.00  0.00          Ve"
line7 = "HETATM  104 APOL STP C   3     117.942 -80.124 -54.380  0.00  0.00          Ve"
# Correct the formatting of the line
print(line_formatting(line))
print(line_formatting(line2))
print(line_formatting(line3))
print(line_formatting(line4))
print(line_formatting(line5))
print(line_formatting(line6))
print(line_formatting(line7))

