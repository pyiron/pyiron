import os

command = list()
command.append("sphinx-apidoc -f -o ./apidoc ../pyiron/")
command.append("make html")


os.chdir(os.getcwd())
os.system(command[0])
os.system(command[1])