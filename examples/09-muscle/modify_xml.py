import re

# Name of your input and output files
input_filename = "muscles.xml"
output_filename = "muscles_modified.xml"

with open(input_filename, 'r') as file:
    xml_data = file.read()

# Function to modify the contents inside <linkName> tags
def update_link_names(match):
    link_name = match.group(1)
    
    # 1. Replace dashes with underscores
    link_name = link_name.replace('-', '_')
    
    # 2. Replace 'hip' with 'pelvis'
    if link_name == 'hip':
        link_name = 'pelvis'
        
    if link_name == 'base':
        link_name = 'pelvis'
        
    if link_name == 'chest':
        link_name = 'torso'
        
    return f"<linkName>{link_name}</linkName>"

# Find all <linkName> tags and apply the modification function
modified_xml = re.sub(r'<linkName>(.*?)</linkName>', update_link_names, xml_data)

# Save the modified XML to a new file
with open(output_filename, 'w') as file:
    file.write(modified_xml)

print(f"Successfully modified and saved to {output_filename}")