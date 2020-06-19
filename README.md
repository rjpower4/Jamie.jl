<p align="center">
  <a href="" rel="noopener">
 <img width=200px height=229px src="assets/icon.svg" alt="Project logo"></a>
</p>

<h3 align="center">Jamie.jl</h3>

<div align="center">

  ![Dev](https://github.com/rjpower4/Jamie.jl/workflows/CI/badge.svg?branch=master)
  [![GitHub Issues](https://img.shields.io/github/issues/rjpower4/Jamie.jl.svg)](https://github.com/rjpower4/Jamie.jl/issues)
  [![GitHub Pull Requests](https://img.shields.io/github/issues-pr/rjpower4/Jamie.jl.svg)](https://github.com/rjpower4/Jamie.jl/pulls)
  [![License](https://img.shields.io/badge/license-MIT-blue.svg)](/LICENSE)
  [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://rjpower4.github.io/Jamie.jl/stable)
  [![Dev](https://img.shields.io/badge/docs-latest-blue.svg)](https://rjpower.github.io/Jamie.jl/latest)

</div>

---

<p align="center"> Simple Astrodynamics Analysis Toolkit Written in Julia
    <br> 
</p>

## 📝 Table of Contents
- [About](#about)
- [Getting Started](#getting_started)
- [Usage](#usage)
- [Built Using](#built_using)
- [Contributing](#contributing)
- [Authors](#authors)

## 🧐 About <a name = "about"></a>
`Jamie` is a library primarily intended for use in scripts and packages for astrodynamics analyses. 
A key goal of `Jamie` is to be *simple* and not burden the user with a complex API.

## 🏁 Getting Started <a name = "getting_started"></a>
and the required dependencies and packages will be installed.

### Prerequisites
This software requires the [Julia Programming Language](https://julialang.org/) which is free and open source software.

```
Give examples
```

### Installing
Getting installing `Jamie` is as simple as installing it through Julia's package management system.
Currently, the package is **not** registered with Julia, so to install simply open up a Julia REPL and type

```
julia> using Pkg

julia> Pkg.add(PackageSpec(url="https://github.com/rjpower4/Jamie.jl.git"))
```

and you're ready to go!
Check it out by grabbing some Lagrange points

```
julia> using Jamie

julia> equilibrium_solutions(CrtbpSystem(0.012))
```

## 🔧 Running the tests <a name = "tests"></a>
To run the tests for `Jamie`, enter the package management mode in the REPL via the `]` key and type

``` 
(@1.4) pkg> test Jamie
```

## 🎈 Usage <a name="usage"></a>
Add notes about how to use the system.

## 🚀Contributing <a name = "contributing"></a>

## ✍️ Authors <a name = "authors"></a>
- [@rjpower4](https://github.com/rjpower4) - Idea & Initial work

See also the list of [contributors](https://github.com/kylelobo/The-Documentation-Compendium/contributors) who participated in this project.
