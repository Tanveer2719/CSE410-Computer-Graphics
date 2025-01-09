# CSE-410: Computer Graphics Sessional

Welcome to the **CSE-410: Computer Graphics Sessional** repository! This repository contains implementations of various offline assignments completed as part of the course. The assignments involve building foundational graphics concepts, implementing key transformations, and developing advanced rendering techniques such as ray tracing.

---

## Table of Contents

1. [Introduction](#introduction)
2. [Assignments](#assignments)
   - [Offline 1: Camera and Object Simulation](#offline-1-camera-and-object-simulation)
   - [Offline 2: Raster-Based Graphics Pipeline](#offline-2-raster-based-graphics-pipeline)
   - [Offline 3: Ray Tracing](#offline-3-ray-tracing)
3. [How to Run](#how-to-run)


---

## Introduction

**CSE-410: Computer Graphics Sessional** is one of the most challenging and exciting sessionals at BUET. The assignments in this course are designed to provide hands-on experience with core computer graphics concepts, including camera control, transformations, rasterization, and ray tracing.

---

## Assignments

### Offline 1: Camera and Object Simulation
Path: [Offline 1](./Offline%201)

In this assignment, the focus was on implementing a fully controllable camera and object interactions. The key features include:

- **Camera Controls**: Implemented movement, rotation, and tilting of the camera.
- **Controllable Rolling Ball**: Simulated the movement of a ball synchronized with the camera's position.
- **Magic Cube**: Created a cube that morphs from an octahedron to a sphere and vice versa, based on the camera's position.

### Offline 2: Raster-Based Graphics Pipeline
Path: [Offline 2](./Offline%202)

This assignment involved developing a **raster-based graphics pipeline** similar to the one used in OpenGL. The pipeline consists of six stages, of which the following four were implemented:

1. **Modeling Transformation**: Transforming objects from their local space to world space.
2. **View Transformation**: Converting world space coordinates to camera space.
3. **Projection Transformation**: Implementing perspective projection to map 3D coordinates onto a 2D screen.
4. **Clipping & Scan Conversion**: Using the Z-buffer algorithm to perform clipping and rendering.

### Offline 3: Ray Tracing
Path: [Offline 3](./Offline%203)

The final assignment focused on **ray tracing**, building upon the implementations of Offline 1 and 2. The features include:

- **Phong Lighting Model**: Implemented realistic illumination with ambient, diffuse, and specular components.
- **Recursive Reflection**: Added support for reflective surfaces using recursive ray tracing.
- **Memory Management**: Efficiently handled large data structures and objects for rendering.

---

## How to Run

To run the assignments:

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/cse-410-computer-graphics-sessional.git
   cd cse-410-computer-graphics-sessional
   ```

2. Navigate to the desired assignment folder (e.g., `Offline 1`):
   ```bash
   cd "Offline 1"
   ```

3. Compile and run the project:
   ```bash
   g++ main.cpp -o output -lGL -lGLU -lglut
   ./output
   ```

4. Follow the on-screen instructions to interact with the graphics application.

---


Feel free to explore the repository and enhance the implementations. If you have any questions or suggestions, don't hesitate to open an issue or contact me!

