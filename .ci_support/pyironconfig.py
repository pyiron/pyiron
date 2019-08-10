import os


def main():
    current_path = os.path.abspath(os.path.curdir)
    top_level_path = current_path.replace('\\', '/')
    resource_path = os.path.join(current_path, "tests", "static").replace('\\', '/')
    if os.name == 'nt':
        top_level_path = top_level_path[0].upper() + top_level_path[1:]
        resource_path = resource_path[0].upper() + resource_path[1:]
    pyiron_config = os.path.expanduser('~/.pyiron').replace('\\', '/')
    if not os.path.exists(pyiron_config):
        with open(pyiron_config, 'w') as f:
            f.writelines(['[DEFAULT]\n',
                          'TOP_LEVEL_DIRS = ' + top_level_path + '\n',
                          'RESOURCE_PATHS = ' + resource_path + '\n'])
    else:
        print('config exists')


if __name__ == '__main__':
    main()