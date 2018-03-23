{
    "name" : "pyiron",
    "sub-menu" : [
        {
            "name" : "import",
            "snippet" : ["import numpy as np",
                         "%matplotlib notebook",
                         "import matplotlib.pylab as plt",
						 "from pyiron.project import Project"
                        ]
        },
        {
            "name" : "pyiron path",
            "snippet" : ["import sys",
                         "sys.path.append('/u/neugebau/PyIron')"]
        },
        {
            "name" : "create project",
            "snippet" : ["pr = Project('your_name')"]
        }
    ]
}