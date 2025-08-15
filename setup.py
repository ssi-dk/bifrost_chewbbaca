from setuptools import setup, find_namespace_packages

setup(
    name="bifrost_chewbbaca",
    version="2.1.0",
    url="https://github.com/ssi-dk/bifrost_chewbbaca",
    author="Kim Ng",
    author_email="kimn@ssi.dk",
    license="MIT",
    packages=find_namespace_packages(),
    python_requires=">=3.6",
    package_data={"bifrost_chewbbaca": ["config.yaml", "pipeline.smk"]},
    include_package_data=True,
    install_requires=["bifrostlib >= 2.0.11"],
)
