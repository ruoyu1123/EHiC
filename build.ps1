$ErrorActionPreference = "Stop"

$compiler = if ($env:CXX) { $env:CXX } else { "g++" }
$output = "hicreate.exe"
$sources = @(
    "main.cpp",
    "matrix.cpp",
    "reference.cpp",
    "fragmenter.cpp",
    "simulator.cpp",
    "long_read_common.cpp",
    "porec.cpp",
    "cifi.cpp"
)

& $compiler "-std=c++17" "-O2" "-Wall" "-Wextra" "-pedantic" "-pthread" "-o" $output @sources
