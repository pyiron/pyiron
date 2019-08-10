import os


def main():
    current_path = os.path.abspath(os.path.curdir)
    top_level_path = current_path.replace('\\', '/')
    resource_path = os.path.join(current_path, "tests", "static").replace('\\', '/')
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