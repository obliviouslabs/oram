-------------------------------
Go Binding Example
-------------------------------
This example shows how to call the C++ codes from Rust.

1. Build the example with $ rm -r target && cargo build
2. Run unit tests with $ cargo test

-------------------------------

For efficiency and obliviousness consideration, the lengths of keys and values in the oblivious map are fixed at compile time. Since the Rust compiler cannot compile C++ code, user needs to declare the desired lengths in omap/interface/omap_declare.cfg. For instance, to define an omap with 20-byte key and 32-byte value, add
DECLARE_OMAP(20, 32)
to omap_declare.cfg, and the omap can be referred to as OMapBinding_20_32 in the rust program.
