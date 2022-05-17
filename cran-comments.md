## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)
- R-hub linux-x86_64-rocker-gcc-san (r-devel)

## R CMD check results
> On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Sebastian Hellmann <sebastian.hellmann@ku.de>'
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    Busemeyer (17:75)
    DSD (17:54, 19:33)
    Hellmann (15:46)
    Henrik (20:8)
    Pleskac (17:65)
    Singmann (20:15)
    al (15:58)
    dynWEV (17:35, 19:21)
    et (15:55)
    preprint (16:9)
    rtdists (19:53)
  
  Found the following (possibly) invalid DOIs:
    DOI: 10.1037/a0019737
      From: DESCRIPTION
      Status: Forbidden
      Message: 403

> On windows-x86_64-devel (r-devel)
  checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

0 errors √ | 0 warnings √ | 2 notes x

## Comments
All words identified as possibly misspelled are checked and are either part of cited authors' names or confidence models. I checked again the DOI and get to the right article page (also from a different device). 

I can't imagine where the note about the file 'lastMiKTeXException' comes from since no such file nor any TeX file is part of my package. 


## Additional rhub checks

The above checks were the results from the `rhub::check_for_cran()` call. I ran additional checks with the source package file:

- rhub::check_for_cran(path="../dynConfiR_0.0.1.tar.gz", platforms="macos-highsierra-release-cran")
  
  (macOS 10.13.6 High Sierra, R-release, CRAN's setup)
  
  Result: no errors, warnings, or notes
- rhub::check_for_cran(path="../dynConfiR_0.0.1.tar.gz", platforms="macos-m1-bigsur-release")
  
  (Apple Silicon (M1), macOS 11.6 Big Sur, R-release)
  
  Result: no errors, warnings, or notes

- The call rhub::check(platform="macos-m1-bigsur-release") resulted in a fail to build the vignette, though. I think this was due to a problem with the ggplot2 package. Still, I could not replicate the error on my own device, where the check did not produce errors (MacPro (from 2013) with BigSur (v. 11.6.3) with R 4.2.0 and ggplot2 (3.3.6))

