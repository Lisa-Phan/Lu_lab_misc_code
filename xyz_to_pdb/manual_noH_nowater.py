link = r"C:\Users\lisap\Lab_code\xyz_to_pdb\MCPB_OUT_solv_noH_nowater.pdb"

new_string = ''
with open(link, 'r') as file: 
    for line in file:
        all_fields = line.split()
        if any(field in all_fields for field in ['WAT', 'H',
                                                 'H1', 'H2', 'H3',
                                                 'HA', 'HA1'
                                                 'HB', 'HB2', 'HB3',
                                                 'HD', 'HD1', 'HD2', 'HD3', 
                                                 'HG', 'HG1',
                                                 'HIP', 'HID', 'HM1']):
            continue 
        else:
            print(line)
            new_string += line

with open(r"C:\Users\lisap\Lab_code\xyz_to_pdb\MCPB_OUT_solv_noH_nowater.pdb", 'w') as file:
    file.write(new_string)

