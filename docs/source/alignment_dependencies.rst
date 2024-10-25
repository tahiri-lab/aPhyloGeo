.. _alignment_dependencies:

Alignment Dependencies Installation
====================================

To use the **MUSCLE**, **CLUSTALW**, and **MAFFT** alignment methods, please ensure the following dependencies are installed. To change the alignment method, modify the `alignment_method` parameter in the `params.yaml` file.

1. **MUSCLE**:
   - Download the appropriate version from the MUSCLE 5.1.0 release:
     - `MUSCLE 5.1.0 Release <https://github.com/rcedgar/muscle/releases/tag/5.1.0>`_
   - **Windows**:
      - Download `muscle5.1.win64.exe <https://github.com/rcedgar/muscle/releases/download/5.1.0/muscle5.1.win64.exe>`_
   - **Linux**:
      - Download `muscle5.1.linux_intel64 <https://github.com/rcedgar/muscle/releases/download/5.1.0/muscle5.1.linux_intel64>`_
   - Place the downloaded file in the `aphylogeo/bin` directory.

2. **CLUSTALW**:
   - Download the latest version from the official website:
     - `CLUSTALW Download <http://www.clustal.org/download/>`_
   - Follow the installation instructions.
   - Ensure the CLUSTALW executable is in the `aphylogeo/bin` directory.

3. **MAFFT**:
   - Download the latest version from the official website:
    - `MAFFT Download <https://mafft.cbrc.jp/alignment/software/windows_without_cygwin.html>`_
   - Follow the installation instructions.
    - **Windows**:
        - Place the **mafft-win** folder in the `aphylogeo/bin` directory.
    - **Linux**:
        - Create a folder named **mafft-linux64** in the `aphylogeo/bin` directory.
        - Move the **mafft.bat** file into the **mafft-linux64** folder.

Ensure that all executables are accessible in the specified directories before running the alignment methods.