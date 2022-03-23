# used for finding the corresponding package

import sasktran as sk
import imp
print(imp.find_module('sasktran'))

# used for finding a global configuration file contains some settings, particularly the location of large databases

print(sk.config.user_config_file_location())


