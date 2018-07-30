from pyiron_base.objects.generic.hdfio import FileHDFio
from pyiron_base.objects.generic.parameters import GenericParameters
from pathlib2 import Path


class Notebook(object):
    """
    class for pyiron notebook objects 
    """
    @staticmethod
    def get_custom_dict():
        folder = Path('.').cwd().parts[-1]
        hdf_file = Path('.').cwd().parents[1]/folder
        hdf_file = str(hdf_file)+'.h5'
        if Path(hdf_file).exists():
            hdf = FileHDFio(hdf_file)
            custom_dict = GenericParameters()
            for k,v in zip(hdf[folder+'/input/custom_dict/data_dict']['Parameter'], hdf[folder+'/input/custom_dict/data_dict']['Value']):
                custom_dict[k] = v
            return custom_dict
        else:
            print(hdf_file, 'not found')
            return None

