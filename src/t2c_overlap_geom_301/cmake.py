import os

files = []

for entry in os.scandir('/Users/rinkevic/Development/Litmus/build'):
    if entry.is_file():
        if entry.name[-4:] == ".cpp":
            files.append(entry.name)

files.sort()

with open('CMakeLists.txt', 'w') as f:
    # write wrapper start
    f.write('target_sources(vlxobjs\n')
    f.write('  PRIVATE\n')
    for file in files:
        f.write(4 * ' ' + file + '\n')
    # write wrapper end
    f.write(')\n\n')
    f.write('target_include_directories(vlxobjs\n')
    f.write(' PUBLIC\n')
    f.write('    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>\n')
    f.write(' )\n')
