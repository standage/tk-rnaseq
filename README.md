## tk-rnaseq:  RNA-Seq analysis toolkit
Daniel S. Standage: <daniel.standage@gmail.com>

**17 Jul 2014**: After relying on the scripts in this toolkit extensively for the last several months, I can gladly say that the time investment required to clean up and document these scripts was well worth it and surely saved me hours upon hours of time had I not been as diligent. That's not to say this level of craftsmanship is worthwhile for *every* piece of code I write, but for anything I know I will use again I never regret the time I've spent improving the interface or the documentation.

Just a note: this toolkit was designed specifically to handle data produced by RSEM and EBSeq. I expect that it shouldn't be too hard to adapt the toolkit for workflows utilizing alternative tools, though. In particular, most R/Bioconductor packages that do the differential expression calls seem to have settled on a *de facto* input format, which is nice.

**14 Nov 2013**: A project I am working on has an RNA-Seq component that has demanded my attention for several months now. In pursuit of this project, I have been copying data back and forth between several machines, leaving a disheveled pile of half-baked READMEs, scripts, and procedures in my wake. I'm OK with this--most of the little scripts I write as part of my research will never again see the light of day. But there are some that I keep re-using and re-adapting over and over again.

This package is my effort to centralize the most important RNA-Seq analysis scripts and procedures that I keep coming back to: the ones worth documenting and generalizing so that I or anyone else can come back in 3-6 months, know what the heck everything is and does, and use them for something else.

Rather than creating a directory and README for each tool, for now I am placing each tool's documentation in that tool's source code. I'm expect this will be a satisfactory solution for most of the tools I store in the package.

