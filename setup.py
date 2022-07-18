from setuptools import find_packages, setup

extra_index_urls = []
packages = []

with open("requirements.txt", encoding="utf-8") as file:
    for line in map(str.strip, file):
        if line:
            if line.startswith("-f"):
                extra_index_urls.append(line.split()[1])
            else:
                packages.append(line)


setup(
    name="sars_processing",
    version="0.0.3",
    packages=find_packages(),
    description="Process sars-cov-2 full genome sequences and get mutations from builded phylogenetic tree",
    author="kpotoh",
    license="MIT",
    install_requires=packages,
    dependency_links=extra_index_urls,
    python_requires=">=3.8",
)
