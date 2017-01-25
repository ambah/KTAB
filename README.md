[![Build Status](https://travis-ci.org/KAPSARC/KTAB.svg?branch=master)](https://travis-ci.org/KAPSARC/KTAB)
[![Build status](https://ci.appveyor.com/api/projects/status/iebb8gg63ugqwfpi/branch/master?svg=true)](https://ci.appveyor.com/project/davidrpugh/ktab/branch/master)


# KTAB #


KTAB is an open-source toolkit for assembling agent-based models of negotiation and bargaining. The main development site for KTAB is the King Abdullah Center for Petroleum Studies and Research, [KAPSARC](http://www.kapsarc.org). Research is distributed at [OpenKAPSARC](https://www.kapsarc.org/openkapsarc/) and the main page for KTAB information is [here](https://www.kapsarc.org/openkapsarc/kapsarc-toolkit-for-behavioral-analysis-ktab/).

All KTAB models represent stochastic decision-making among comparatively small numbers of stakeholder groups (roughly 5 to 50), within the paradigm of "Probabilistic Condorcet Elections". PCEs are used to estimate the likelihood of different outcomes from a collective decision making process, depending on what kinds of options each agent has or can generate, how they value those options, and what kinds of coalitions they can form to support or oppose each option. Details can be found in the online documentation as well as the main KTAB  [page](https://www.kapsarc.org/openkapsarc/kapsarc-toolkit-for-behavioral-analysis-ktab/).

## Examples ##

Examples are provided to illustrate the range of modeling which can be done in the framework and to provide templates for further work. Examples in the kmodel directory include the following:

- Using an input/output economic model to bargain over tax/subsidy rates,
- Simple discrete matching
- Negotiating the priority order to tasks to undertake (and which to drop)


- A simple but highly parameterizable version of the [spatial model of politics](examples/smp/README.md). The command line version is `smpc`;  GUI versions will be called `smpg` and `smpq`. Dummy data is provided to illustrate the format for CSV input files.


- A simple model of bargaining over the order in which projects should be addressed, taking into account budget limits. Dummy data is provided.


- A simple example of using KTAB to fit parameters. Similar to econometric models, the parameters of a PCE are adjusted until the result closely matches the expected results. The fitted parameters could then be used to analyze possible changes.

These are merely illustrative examples.  For rigorous work, one could make whatever local modifications were needed to produce a model which met one's own standards for validation.

## Building KTAB ##

KTAB is written in portable, cross-platform C++11.

KTAB uses CMake for configuration; downloads and instructions can be obtained from [www.cmake.org](http://www.cmake.org). The basic procedure is to configure kutils with CMake, then build the library and examples with your favorite IDE. Then configure kmodel with CMake and build it. Detailed build instructions can be found in the subdirectories, starting with kutils' [README](KTAB/kutils/README.md). After building kutils, you can follow kmodel's [README](KTAB/kmodel/README.md). Following these instructions, under both Windows and Linux, it is not necessary to modify your PATH at any point of the installation.

Example of how to build models outside the main library structure, as one would do for a custom model, are in the examples directories; see the overall [README](examples/README.md) and the README files for each.


The code is compiled, linked, and tested with absolutely no modification, simultaneously on both Windows and Linux: there is one single body of cross-platform code, not  separate versions for each OS. All platform- or system-specific configuration is done by CMake. Because no platform-specific code is allowed, porting to Mac or other platforms should be feasible. It was developed and tested on 64-bit Windows and Linux platforms, including the following:

* Windows
  * Windows 7 Professional with Visual Studio 2010 Express,
  * Windows 7 Professional with Visual Studio 2013, Update 5,
  * Windows 8 &amp; 10 Professional with Visual Studio 2015 Community,
  * Windows Server 2012 with Visual Studio 2015 Community,
* Linux
  * Fedora 22 with gcc/g++ 5.1,
  * openSUSE 13.2 with clang/clang++ 3.5,
  * Debian 8.1 with gcc/g++ 4.9,
  * Ubuntu 16.04 with gcc/g++ 5.4.

### Language Versions ###

The C++11 standard is supported by the following versions of
these common compilers:

* GCC 4.9.0, or higher
* Clang 3.5.0, or higher
* Visual Studio 2013, or higher

We expect to upgrade to C++17 when it becomes well-supported
by gcc, clang, and VS, especially structured bindings, aka
parallel assignment with pattern matching.

## User Interfaces ##

A graphical front-end is under development. The current GUI version of the SMP, smpq, has been developed in [Qt](https://www.qt.io/). Qt is dual-licensed under both commercial and open-source terms; details can be found [here](https://www.qt.io/qt-licensing-terms/). The Qt licensing page states that dynamically linking an application to the open-source Qt libraries does not make the application open-source.

Previous work toward a graphical front-end was based on the FLTK project (http://www.fltk.org) [FLTK](http://www.fltk.org). FLTK is a cross-platform, C++ tookit provided under the terms of the GNU Library Public License, Version 2 with exceptions that allow for static linking; details can be found [here](http://www.fltk.org/COPYING.php).  The source code and build instructions for the most recent stable version can be obtained  [here](http://www.fltk.org/software.php). FLTK also uses CMake for configuration, then your favorite IDE to build.

Because the libraries that analyze negotiation are completely separate from the applications that use them, developers are perfectly free to write command line applications or use other GUI toolkits. Along these lines, the GUI work is expected to diversify into other toolkits, such as [WxWidgets](https://www.wxwidgets.org/), or even a web interface such as [Wt](http://www.webtoolkit.eu/wt).

For authoritative ruling on licensing issues, consult appropriate counsel.

## Data Input Formats ##

To read and write XML data files, KTAB uses the [TinyXML2](https://github.com/leethomason/tinyxml2) library, which is available on [GitHub](https://github.com/).


## Contributing and License Information ##


If you are interested in contributing code, ideas, or
data to KTAB, please contact ktab@kapsarc.org


KTAB is released under The MIT License (Expat).
For details, see the following URLs:

- [http://opensource.org/](http://opensource.org/)
- [http://opensource.org/licenses/MIT](http://opensource.org/licenses/MIT)


----------

Copyright KAPSARC. Open source MIT License.

----------

<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-51793176-2', 'auto');
  ga('send', 'pageview');

</script>
