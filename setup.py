import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
    setuptools.setup(
        name='sequtils',
        version='0.1',
        scripts=[],
        author="Jan Christian Refsgaard",
        author_email="jancrefsgaard@gmail.com",
        description="Biological Sequence Utility Scripts",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/jancr/sequtils",
        #  packages=setuptools.find_packages(),
        packages=['sequtils'],
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
    )
