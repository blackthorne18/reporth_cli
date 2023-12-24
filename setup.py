from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = fh.read()

setup(
    name='repinclusterer',
    version='0.7.0',
    author='Prajwal Bharadwaj',
    description='Cluster DNA Sequences based on orthology',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/blackthorne18/repinclusterer_cli',
    py_modules=['rc_entry', 'mkclus'],
    packages=find_packages(),
    install_requires=[requirements],
    python_requires='>=3.10',
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],
    entry_points='''
        [console_scripts]
        repinclusterer=rc_entry:main
    '''
)
