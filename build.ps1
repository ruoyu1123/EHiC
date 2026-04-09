$ErrorActionPreference = "Stop"

$compiler = if ($env:CXX) { $env:CXX } else { "g++" }
$output = "hicreate.exe"
$sources = @(
    "main.cpp",
    "matrix.cpp",
    "reference.cpp",
    "fragmenter.cpp",
    "simulator.cpp"
)

& $compiler "-std=c++17" "-O2" "-Wall" "-Wextra" "-pedantic" "-o" $output @sources
