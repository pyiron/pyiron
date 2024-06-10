import matplotlib.pyplot as plt
import numpy as np
from pyiron import Project


if __name__ == "__main__":
    pr = Project('thermo')
    basis = pr.create.structure.bulk(name='Fe', cubic=True, a=2.75)
    line_strain_list = np.linspace(0.95, 1.05, 7)**(1 / 3) - 1
    for strain in line_strain_list:
        dft = pr.create.job.Gpaw(job_name=('sample', strain))
        dft.set_encut(320.0) # Optional, among other plane wave dft parameters
        dft.structure = basis.apply_strain(strain, return_box=True)
        dft.run()
    results = {'energy': [], 'volume': []}
    for strain in line_strain_list:
        dft = pr.load(('sample', strain))
        results['volume'].append(dft.structure.get_volume())
        results['energy'].append(dft.output.energy_pot[-1])
    coeff = np.polyfit(results['volume'], results['energy'], 3)
    equi_volume = np.roots(np.polyder(coeff)).min()
    print('Equilibrium volume:', np.round(equi_volume, decimals=3), 'A^3')
    equi_bulk_mod = np.polyval(np.polyder(coeff, 2), equi_volume) * equi_volume
    print('Equilibrium bulk modulus:', np.round(equi_bulk_mod, decimals=3), 'eV/A^3')
    dft = pr.create.job.Gpaw(job_name='gpaw')
    dft.set_encut(320) # Again, optional
    dft.structure = basis.copy()
    murn = dft.create_job(job_type=pr.job_type.Murnaghan, job_name='murn')
    murn.input['num_points'] = 7
    murn.input['vol_range'] = 0.05
    murn.run()
    fit_dict = {
        "Vinet": murn.fit_vinet(),
        "Murnaghan": murn.fit_murnaghan(),
        "Birch-Murnaghan": murn.fit_birch_murnaghan(),
        "Polynomial": murn.fit_polynomial(),
    }
    print(fit_dict)
