# Package

version       = "0.1.0"
author        = "jcalendo"
description   = "Fast transcript integrity number calculation"
license       = "MIT"
srcDir        = "src"
bin           = @["tin"]


# Dependencies

requires "nim >= 2.0.0"
requires "hts >= 0.3.30"
requires "lapper >= 0.1.8"
requires "cligen >= 1.9.6"