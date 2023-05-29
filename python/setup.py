import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="libtetrabz",
    version="0.1.0",
    author="Mitsuaki Kawamura",
    author_email="kawamitsuaki@gmail.com",
    description="Optimized tetrahedron method to perform brillouin-zone integral",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://libtetrabz.osdn.jp/",
    project_urls={
        "Bug Tracker": "https://osdn.net/projects/libtetrabz/forums/",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    ext_modules=[setuptools.Extension('libtetrabz.libtetrabzc', ['src/libtetrabz/libtetrabz.c'])]
)
