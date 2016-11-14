WCS-ADDPV
=========

wcs-addpv.py is a script to derive PV polynomial WCS correction terms
from SIP corrections and add the PV keywords to a FITS image header,
so that the WCS can be (better) understood by tools in the
SEtractor/Scamp/Swarp suite of software tools.

This is most useful where you have derived WCS corrections in another
way (e.g., through astrometry.net) and the corrections are presented
in the SIP (Simple Imaging Polynomial) format (note the CTYPEi
keywords: these should then be ``RA---TAN-SIP`` and ``DEC--TAN-SIP``).


Usage
=====

In its simplest form, run it as::

        ./wcs-addpv.py somefile.fits

Running it without input files or with the ``--help`` option gives a
quick summary of the script. wcs-addpv can take one or more input FITS
images.

wcs-addpv will alter your input FITS files *in place*, so make a
backup of the file if needed. Note that wcs-addpv only *adds
redundant* keywords: it will not replace or remove keywords, except
existing ``PVi_j`` ones. In particular, it does not change the CTYPEi
keywords. (For the how and why, see the Note in the Background.)

wcs-addpv only uses the primary image, that is, HDU 0. All other HDUs
in the file are left untouched. Thus, it works best on simple single
image FITS files.

Options
-------

Some options require some understanding of the algorithm underlying
the script. See the background section for that.

- ``--order``: the polynomial order for the PV corrections. The
  default is 3.

  The zeroth and first order are never touched, and always left at
  their defaults of 0 and 1, respectively. So the minimum usable
  order is 2, and the default is set to 3.

  A good starting point, in the case of existing SIP corrections,
  is to use the same order as that of the SIP corrections, or
  lower. A higher order is much more likely to lead to erroneous
  corrections.

- ``--radial``: use the radial PV corrections. The default is off.

  This switches on the radial term in the PV corrections. The
  Sextractor/Scamp/Swarp suite ignore the radial terms, so it's
  generally fine to leave this off.

- ``--ndata``: number of data points to use in the fit. Default 64.

  This sets the approximate number of data points to use in the
  fit. The option takes one or two arguments. The latter sets the
  values independently for each dimension x and y, the former uses
  the same value for each dimension.

  The total number of data points is then ndata x ndata. This is
  the approximate amount of data points used in the fit (see the
  background section, the algorithm subsection).

  The number is approximate, as it tries to subdivide the image
  data so that there are ndata points for the respective
  dimension. The default is therefore 64, since most image
  dimension are some multiple of 64.

- ``--verbose``: show a bit more output. Can be used multiple times.
  The default is no information at all.

  Use once to obtain output information on the improvement in the
  PV correction before and after the fit. Use twice to obtain
  information at the debugging level.



Installation
============

Installation is simply::

    python -m pip install git+https://github.com/evertrol/sippv.git


You can also grab the standalone ``wcs-addpv.py`` script from the
``scripts`` subdirectory. This standalone script will not be installed
when installing with ``pip``; instead, ``pip`` will install a
``sippv-add`` script that wraps around the module.


Requirements
------------

wcs-addpv requires the following libraries:

- numpy
- scipy
- astropy

wcs-addpv has only been tested with the current astropy version (1.2),
but I expect it to run with any recent 1.x astropy wcs. Similar for
numpy and scipy: any recent version should work.

It has been used mostly with Python 3.5, but has been tested with 2.7
as well.

Alternatives
============

As an alternative to fixing the headers, there are a few alternatives.

- For SExtractor: convert the pixel coordinates in the catalog file to
  proper world coordinates::

    with astropy.io.fits.open("somefile.fits") as hdulist:
        header = hdulist[0].header
        wcs = astropy.wcs.WCS(header)

    table = astropy.table.Table.read("sex.cat", format="ascii.sextractor")
    xpix = table['X_IMAGE']
    ypix = table['Y_IMAGE']
    pixcoords = np.asarray([xpix, ypix]).T
    # all_pix2world() includes the SIP corrections
    wcoords = wcs.all_pix2world(pixcoords, 1).T
    table['ra'] = wcoords[0]
    table['dec'] = wcoords[1]

- A simple (not very well tested) alternative to ``SWarp`` is to use
  the utilities in ``scipy.ndimage``::

    with astropy.io.fits.open("somefile.fits") as hdulist:
        inhead = hdulist[0].header
        inwcs = astropy.wcs.WCS(inhead)
        inshape = hdulist[0].data.shape
    with astropy.io.fits.open("reffile.fits") as hdulist:
        refhead = hdulist[0].header
        refwcs = WCS(refhead)
        refdata = hdulist[0].data.copy()

    xpix, ypix = np.mgrid[0:inshape[1], 0:inshape[0]]
    coords = np.asarray([xpix.ravel(), ypix.ravel()]).T
    wcoords = inwcs.all_pix2world(coords, 1, ra_dec_order=True)
    refcoords = refwcs.all_world2pix(wcoords, 1, ra_dec_order=True)

    refdata = refdata.T
    output = map_coordinates(refdata, refcoords.T)
    output = output.reshape(refdata.shape).T

    # The output file should have the same WCS as the input file
    # relax=True keeps the RA/DEC--TAN-SIP and polynomial coeffs keys
    outhead = inhead.copy()
    outhead.update(refwcs.to_header(relax=True))
    fits.PrimaryHDU(data=output, header=outhead).writeto(
        "newfile.fits", clobber=True)



Background
==========

In (difference) imaging pipelines, where astrometry, photometry and
image subtraction is used, commonly used tools are astrometry.net,
SExtractor, SWarp and HOTPANTS. The SExtractor and SWarp tools rely on
the WCS to obtain world coordinates. Unfortunately, they use a
different way (TPV) to interpret distortions to the WCS than is
nowadays more commonly used (SIP; used ine.g., in astrometry.net). As
such, the WCS distortions calculated by astrometry.net are not applied
to the output of SExtractor and SWarp, leading to incorrect positions
and image alignments.

History
-------

The undistorted WCS system is well described and an accepted standard.
For distortions, no standard currently exists (to my knowledge), and
historically, two alternatives were available. One is the TPV system,
where polynomial distortions are added to the intermediate world
coordinates.

The intermediate world coordinates (Calabretta & Greisen, 2002 [#f1]_ ) are
given by the initial offset and rotation/shear transformation::

    | x1 | = | cd1_1  cd1_2 |  | p1 - crpix1 |
    | x2 | = | cd1_2  cd2_2 |  | p2 - crpix2 |


When no distortions are present, these intermediate world coordinates
are transformed to final world coordinates (right ascension and
declination), according to their projection type. See Calabretta &
Greisen (2002), section 7.3 for some examples of this.

When distortions are present (both in the image and through correction
terms in the WCS), polynomials corrections are added to the `xi` world
coordinates, of the form::

    x1' = f1(x1, x2, r)
    x2' = f2(x1, x2, r)

where ``r`` is a shorthand for ``r = sqrt(x1*x1 + x2*x2)``, and is set to
0 in the SExtractor/SWarp case. ``f1`` and ``f2`` are essentially
polynomials, similar to the SIP polynomials described below. The
polynomial terms ``Pi_j`` are then stored in a FITS header as ``PVi_j``
keys.

With the corrected intermediate world coordinates, the usual
transformation to actual RA and Dec is then applied.

This applies only to the ``TAN`` (gnomic) projection, and the ``CTYPi``
keywords should be set to ``TPV`` when using these corrections.


The full description of the terms is availabe at
http://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html .

Since this suggested correction type was available early on,
SExtractor and SWarp have incorporated it.

Over time, however, another correction type has gradually become the
de facto standard: SIP, the Simply Imaging Polynomial, corrections.

These corrections are very similar to the PV correction, albeit
without a radial term, but are applied *before* the initial
transformation to intermediate world coordinates::

    | x1 | = | cd1_1  cd1_2 |  | p1 + f1(p1, p2) |
    | x2 | = | cd1_2  cd2_2 |  | p2 + f2(p1, p2) |

with the correction functions polynomials::

    fi(p1, p2) = \Sigma_j_k P_i_j_k p1^j p2^k

When using SIP correction, the ctype keyword is extended to include
``-SIP`` at the end.

More details at http://fits.gsfc.nasa.gov/registry/sip/shupeADASS.pdf



Correction conversion
---------------------

It is possible to convert one type of corrections to the other, albeit
not that straightforward. In fact, this is something the PTF (Palomar
Transient Factory) group has done, as they deal with large images
(thus causing distortions) and a pipeline using Astrometry.net and
SExtractor (Shupe et al, 2012 [#f2]_).

Unfortunately, while they describe the method in their paper, they
have only implemented this in code for a fixed correction order. They
have made their software available [#f3]_, [#f4]_, [#f5]_, but only in
binary format, as a standalone executable.

This situation felt rather unsatisfactory to me. Hence I decided to
implement a more generic solution to this problem, albeit not as
direct as the PTF solution.


Algorithm
---------

The algorithm in wcs-addpv relies on minimizing the difference in the
*known* SIP-corrected coordinates, and the coordinates calculated fro
the to-be-estimated PV terms. Thus, it becomes somewhat of a forward
fitting process.

``astropy.wcs.WCS`` will read the available WCS with the SIP
corrections. It will then obtain the world coordinates for a subset of
image pixels (centers), spread evenly across the image (it is possible
to use all pixels, though this would increase the calculation time.
Currently, about 5000 pixels are used; this is the ``ndata`` parameter
squared).

``wcs-addpv`` then calculates the ``PV`` distortions, starting with a
default set of polynomial coefficients (all zero, except for ``PV1_1``
and ``PV2_1``), and up to a given order (set by the ``order`` keyword).
Note that the zeroth order term is always kept at zero, and the first
order term is always kept at one, since these can essentially be
incorporated in the ``CRVALi`` values and transformation matrix. The
radial correction is normally fixed at zero, but can be fit for as
well, when the ``radial`` option is ``True``.

The sum of the square of the differences in RA and Dec between the two
solutions is then minimized. This yields the ``PV`` correction terms,
which are then written to the header.


.. Note::

The ``PVi_j`` keywords are added, and the ``CTYPEi`` is not altered. This
relies on the fact that SExtractor and SWarp do not look for a ``PTV``
ctype, but instead rely on the ``PVi_j`` keywords to be present to apply
a correction. This is why the script is called ``wcs-addpv``, and not
``sip2pv`` or similar.

This means that a FITS header can have valid SIP corrections, and
(invalid) PV corrections at the same time.

When reading such a header with ``astropy.wcs.WCS``, a warning is issued
and the PV keys are removed from the WCS (since they are redundant,
and have caused issues in the past [#f6]_ ). wcs-addpv silences this
specific warning for convenience, since the PV keywords will be are
intended to be replaced.


Current implementation
----------------------

The current implementation is straightforward, simple but slow. It can
take several seconds to a minute or more to minimize the difference,
depending on the amount of pixels used in the calculation. Lowering
the amount of pixels (``ndata``) will speed up the calculation, but may
reduce the accuracy of the fit.

One of the future plans is to speed up the PV correction section in
the code using either C or Cython.


License
=======

``sippv`` is distributed under the ISC license.

``sippv`` is copyright `the Gravitational-wave Optical Transient
Observer (GOTO) Observatory <http://goto-observatory.org/>`_.


References & footnotes
======================

.. [#f1] http://adsabs.harvard.edu/abs/2002A%26A...395.1077C

.. [#f2] http://web.ipac.caltech.edu/staff/shupe/reprints/SIP_to_PV_SPIE2012.pdf

.. [#f3] astrometry mailing list thread:
         https://groups.google.com/forum/#!topic/astrometry/qrNW5KtsNCc

.. [#f4] ``sip2pv`` and ``pv2sip`` binaries at
         http://data.astrometry.net/pv2sip-binaries-intel/ . These are
         pre-compiled, so will not work on every platform, and require
         ``solve-field`` to be run with ``--tweak-order 4``.

.. [#f5] Astrometry.net ships with a ``wcs-pv2sip`` utility, which is
         probably what is used with the ``--scamp`` option to
         ``solve-field``. I have not yet inspected this further, but
         the naming suggests this utility works the inverse way of
         what I want .

.. [#f6] https://github.com/astropy/astropy/issues/299
