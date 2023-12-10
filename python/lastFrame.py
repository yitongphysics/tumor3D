def keep_last_block(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Find the start of the last block
    last_block_start = None
    for i in range(len(lines) - 1, -1, -1):
        if lines[i].strip() == "NEWFR":
            last_block_start = i
            break

    # If a block was found, write only the last block back to the file
    if last_block_start is not None:
        with open(file_path, 'w') as file:
            file.writelines(lines[last_block_start:])

# Example usage
# keep_last_block('path_to_your_file.txt')

path_file = '/Users/yitongzheng/Documents/Corey/tumor3D/P.pos'
keep_last_block(path_file)

'''
path_file = '/gpfs/gibbs/project/ohern/yz974/tumor3D/invasion/0425_long/test/'

for ii in range(16):
    fstr = path_file + str(ii) + '.pos'
    keep_last_block(fstr)
'''