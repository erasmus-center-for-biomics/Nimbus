Structural variant calling with Nimbus
==========================

Introduction
------------

For structural variant (SV) detection, amplicons are designed over breakpoints. With some adaptations, Nimbus can detect and quantify these SVs. In this guide, the adaptions to the protocol are discussed.

Amplicon design considerations
------------------------------

To effectively call SVs with amplicon based technologies, amplicons must be designed over the break-points. In the case of a typical PCR, primers should be placed on next to both sites of the break point (Figure 1). Furthermore, primers can also be designed to identify the absence of a structural variants. With a clear amplicon design,the presence and absence of amplicons show the presence of structural variants. The analysis of these design can be facilitated by Nimbus.

```amplicon
break-point = //

primer a                  -->
chromosome a: +++++++++++++++++++//++++++++++++++++++++++
primer a*:                            <--

primer b*                   -->
chromosome b: ===================//======================
primer b                              <--

Yields:

if SV:
          +++++++//======
if no SV:
          +++++++//++++++
          =====//======
if "heterozygous" SV:
          +++++++//======
          +++++++//++++++
          =====//======
```

Figure 1 amplicon with a break point
