-------------------------------
Go Binding Example
-------------------------------
This example shows how to compile the C++ codes into a shared library, and called by a Go program. Wrapper codes are generated using the Swig tool.

1. Make sure Go and Swig are installed. Change go.mod with the correct Go version.
2. Switch to the odsl folder.
3. Run $ source gen_binding.sh   (you may need to run with sudo)
4. Run $ go test