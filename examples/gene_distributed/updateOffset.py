from configparser import SafeConfigParser
import collections
# read parameter file
parser = SafeConfigParser()
parser.read('ctparam')

config = collections.namedtuple('Config', 'checkpointFrequency ncombi')

config.ncombi = int( parser.get('ct','ncombi') )

with open('ginstance/offset.txt', 'r') as f:
    offset = int(next(f).split()[1])
with open('ginstance/offset.txt', 'w') as f:
    new_offset_string = "offset " + str(int(offset + config.ncombi)) + "\n"
    f.write(new_offset_string) 

