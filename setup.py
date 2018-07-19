import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='rnafbinv', version='2.0b1',
                 description="Fragment based RNA designer",
                 url="https://github.com/matandro/RNAsfbinv",
                 long_description=long_description,
                 long_description_content_type="text/markdown",
                 author="Matan Drory Retwitzer",
                 author_email="matandro@post.bgu.ac.il",
                 packages=setuptools.find_packages(),
                 classifiers=(
                     "Programming Language :: Python :: 3",
                     "License :: OSI Approved :: MIT License",
                     "Operating System :: OS Independent",
                 )
                 )
