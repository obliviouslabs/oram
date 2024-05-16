%module odsl
%{
#include "stdint.h"
#include "interface/common_impl.hpp"
#include "interface/omap_generic_impl.hpp"
#include "interface/par_omap_impl.hpp"
#include "interface/recoram_impl.hpp"
%}

/* Let's just grab the original header file here */
%include <stdint.i>
%include "interface/common_interface.hpp"
%include "interface/omap_generic_interface.hpp"
%include "interface/par_omap_interface.hpp"
%include "interface/recoram_interface.hpp"
%template(OMapBinding) OMapGenericBinding<8, 8>;
%insert(cgo_comment_typedefs) %{
#cgo LDFLAGS: -lodsl -lbearssl -lstdc++
%}