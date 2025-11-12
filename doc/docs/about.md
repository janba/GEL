# About PyGEL3D

## What is GEL?

GEL is a C++ library of geometry processing tools primarily intended for computer graphics applications. It provides a comprehensive set of data structures and algorithms for working with polygonal meshes, spatial graphs, voxel grids, and other geometric representations.

## History and Development

GEL has been developed over many years at the Technical University of Denmark (DTU) by Jakob Andreas Bærentzen and contributors. It has been used in research, teaching, and practical applications in computer graphics and geometry processing.

## PyGEL3D: Python Bindings

PyGEL3D provides Python bindings for a subset of GEL's functionality, making geometry processing accessible to Python developers and researchers. It wraps the most commonly used features of GEL while maintaining the performance benefits of the underlying C++ implementation.

### Why PyGEL3D?

- **Ease of Use**: Python interface for complex geometry operations
- **Performance**: C++ backend for computational efficiency
- **Jupyter Integration**: Interactive development and teaching
- **Cross-Platform**: Works on Windows, macOS, and Linux
- **Well-Tested**: Years of development and practical use

## Key Features

### Halfedge Mesh Data Structure
A mature and efficient implementation that enables:
- Constant-time topological queries
- Flexible mesh traversal
- Support for arbitrary polygons
- Rich editing operations

### Advanced Algorithms
- **Local Separator Skeletonization**: Graph-based curve skeleton extraction
- **Inverse Skeletonization**: Face Extrusion Quad (FEQ) mesh generation from graphs
- **Rotation System Reconstruction**: Combinatorial point cloud reconstruction
- **Garland-Heckbert Simplification**: Quality mesh decimation
- **Various Smoothing Methods**: Laplacian, Taubin, anisotropic, and more

### Visualization Tools
- OpenGL-based interactive viewer
- Plotly integration for Jupyter notebooks
- Multiple rendering modes
- Export-friendly (HTML with embedded 3D)

## Use Cases

### Research
- Algorithm prototyping
- Comparative studies
- Result visualization
- Data processing pipelines

### Teaching
- Geometry processing courses
- Computer graphics education
- Interactive demonstrations
- Assignment submission (via Jupyter notebooks)

### Development
- Rapid prototyping
- Pre-processing pipelines
- Mesh analysis tools
- Integration with other Python libraries

## License

GEL and PyGEL3D are available under the MIT License. See the LICENSE.md file in the repository for full details.

## Citation

If you use PyGEL3D in your research, please cite:

```
@software{gel2024,
  author = {Bærentzen, Jakob Andreas},
  title = {GEL: Geometry and Graphics Library},
  year = {2024},
  url = {https://github.com/janba/GEL}
}
```

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

See the CONTRIBUTORS.md file for more information.

## Support

- **GitHub Issues**: Report bugs and request features
- **Documentation**: This documentation and inline code comments
- **Examples**: Check the examples directory in the repository

## Related Projects

- **GEL C++ Library**: The underlying C++ library
- **scipy.spatial**: Alternative spatial data structures in Python
- **trimesh**: Another Python library for mesh processing
- **PyMesh**: Mesh processing library with different focus
- **Open3D**: 3D data processing library

## Acknowledgments

PyGEL3D builds on years of development and contributions from:
- Jakob Andreas Bærentzen (primary developer)
- Students and researchers at DTU
- The broader geometry processing community

Special thanks to all contributors who have helped improve GEL and PyGEL3D over the years.

## Links

- [GitHub Repository](https://github.com/janba/GEL)
- [PyPI Package](https://pypi.org/project/PyGEL3D/)
- [DTU Homepage](https://www.dtu.dk/)
