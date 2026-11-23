from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [
        line.strip()
        for line in fh
        if line.strip() and not line.lstrip().startswith("#")
    ]

setup(
    name='reporth',
    version='2.6.2',
    author='Prajwal Bharadwaj',
    description='Cluster DNA Sequences based on orthology',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/blackthorne18/reporth_cli',
    py_modules=['rc_entry'],
    packages=find_packages(),
    install_requires=requirements,
    python_requires='>=3.10',
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],
    entry_points='''
        [console_scripts]
        reporth=rc_entry:main
    '''
)
