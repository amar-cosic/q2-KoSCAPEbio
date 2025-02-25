from setuptools import setup, find_packages

setup(
    name='koscapebio',
    version="1.0.0",
    packages=find_packages(),
    author="Amar Cosic",
    author_email="amar.cosic995@gmail.com",
    description=" QIIME 2 plugin designed for the exploration and analysis of the Klebsiella oxytoca species complex (KoSC) presence within microbial communities",
    license='BSD-3-Clause',
    url="https://github.com/amar-cosic/q2-KoSCAPEbio",
    entry_points={
        'qiime2.plugins': ['q2_koscapebio=q2_koscapebio.plugin_setup:plugin'],
        'console_scripts': [
            'database_curation=q2_koscapebio.database_curation.database_curation:main',
        ]
    },
    package_data={
        'q2_koscapebio': ['data/*', 'database_curation/bla_database/*'],
    },
    classifiers=[
        "Framework :: Qiime2",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)

