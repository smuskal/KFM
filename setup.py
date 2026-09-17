from setuptools import setup, find_packages

setup(
    name="kfm",
    version="2.0.0",
    description="Command line tools for the Kinase Foundation Model v2 models",
    packages=find_packages(include=["kfm", "kfm.*"]),
    python_requires=">=3.10",
    # Pinned, not ranged: the forests are pickled object graphs and a different
    # scikit-learn either refuses them or scores differently.
    # rdkit is a range, matching requirements.txt: an exact pin is uninstallable
    # on Intel macOS, and the ceiling holds installs to the releases measured to
    # compute identical ligand features.
    install_requires=["scikit-learn==1.7.2", "numpy==2.2.6",
                      "joblib==1.5.3", "rdkit>=2024.3.6,<=2025.9.5"],
    extras_require={"sequences": ["torch==2.9.1", "transformers==4.57.3"]},
    entry_points={"console_scripts": ["kfm = kfm.cli:main"]},
)
