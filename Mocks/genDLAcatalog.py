import os
import glob
import sys
import scipy as sp
import argparse
import fitsio


def DLAcat_quickq(INPUT,OUTPUT):
	# Check if output folder exist, else create it.
	if os.path.exists(OUTPUT):
		print('DLA catalog already exists')

		return
		
	else:
		print('Creating catalog, reading files in:',INPUT)
		fs = glob.glob(INPUT+'spectra-16/*/*/truth*.fits')
		fs2 = glob.glob(INPUT+'spectra-16/*/*/spectra*.fits')
		fs = sp.sort(sp.array(fs))
		fs2 = sp.sort(sp.array(fs2))
		f_radla = sp.array([],dtype=float)
		f_decdla = sp.array([],dtype=float)
		f_zdla = sp.array([],dtype=float)
		f_mockiddla = sp.array([],dtype=sp.int64)
		f_zqsodla = sp.array([],dtype=float)
		f_nhidla = sp.array([],dtype=float)

		for j,f in enumerate(fs):
			h = fitsio.FITS(f)
			try:
				mockiddla=h['DLA_META']['QSOID'].read()
			except:
				mockiddla=[]
			if(len(mockiddla)>0):
				h2 = fitsio.FITS(fs2[j])
				ra = h2['FIBERMAP']['RA_TARGET'][:]
				dec = h2['FIBERMAP']['DEC_TARGET'][:]
				mockid = h2['FIBERMAP']['TARGETID'][:]
				zqso = h['TRUTH']['TRUEZ'][:]
				w = sp.in1d(mockid,mockiddla)
				ra = ra[w]
				dec = dec[w]
				zqso = zqso[w]
				mockid = mockid[w]
				frommockidtora = {mockid[i]:ra[i] for i in range(ra.size)}
				frommockidtodec = {mockid[i]:dec[i] for i in range(ra.size)}
				frommockidtozqso = {mockid[i]:zqso[i] for i in range(ra.size)}
				radla = [ frommockidtora[i] for i in mockiddla ]
				decdla = [ frommockidtodec[i] for i in mockiddla ]
				zqsodla = [ frommockidtozqso[i] for i in mockiddla ]

				f_radla = sp.append(f_radla,radla)
				f_decdla = sp.append(f_decdla,decdla)
				f_zqsodla = sp.append(f_zqsodla,zqsodla)
				f_mockiddla = sp.append(f_mockiddla,mockiddla)

				f_zdla = sp.append(f_zdla,h['DLA_META']['ZDLA'][:])
				f_nhidla = sp.append(f_nhidla,h['DLA_META']['NHI'][:])

				h.close()
				h2.close()

			### Sort
			w = sp.argsort(f_mockiddla)
			for el in [f_radla,f_decdla,f_mockiddla,f_zdla,f_mockiddla,f_mockiddla,f_mockiddla,f_nhidla,f_zqsodla]:
    				el = el[w]

			### Save
		out = fitsio.FITS(OUTPUT,'rw',clobber=True)
		cols=[f_radla,f_decdla,f_mockiddla,f_zdla,f_mockiddla,f_mockiddla,f_mockiddla,f_nhidla,f_zqsodla]
		names = ['RA','DEC','DLAID','ZDLA','PLATE','MJD','FIBERID','NHI','ZQSO']
		out.write(cols,names=names)
		out.close()

		
		print('Done creating DLA catalog')

	return 	

def DLAcat_trans(INPUT,OUTPUT):
	# Check if catalog exist
	if os.path.exists(OUTPUT):
		print('DLA catalog already exists')

		return 
		
	else:
		print('Creating catalog. Reading files in:',INPUT)
		fs = glob.glob(INPUT+'/*/*/transmission-16-*.fits')
		fs = sp.sort(sp.array(fs))

		f_radla = sp.array([],dtype=float)
		f_decdla = sp.array([],dtype=float)
		f_zdla = sp.array([],dtype=float)
		f_mockiddla = sp.array([],dtype=sp.int64)
		f_zqsodla = sp.array([],dtype=float)
		f_nhidla = sp.array([],dtype=float)

		for j,f in enumerate(fs):
			h = fitsio.FITS(f)
			ra = h['METADATA']['RA'][:]
			dec = h['METADATA']['DEC'][:]
			zqso = h['METADATA']['Z'][:]
			mockid = h['METADATA']['MOCKID'][:]
			mockiddla = h['DLA']['MOCKID'][:]
			w = sp.in1d(mockid,mockiddla)
			ra = ra[w]
			dec = dec[w]
			zqso = zqso[w]
			mockid = mockid[w]
			frommockidtora = {mockid[i]:ra[i] for i in range(ra.size)}
			frommockidtodec = {mockid[i]:dec[i] for i in range(ra.size)}
			frommockidtozqso = {mockid[i]:zqso[i] for i in range(ra.size)}
			radla = [ frommockidtora[i] for i in mockiddla ]
			decdla = [ frommockidtodec[i] for i in mockiddla]
			zqsodla = [ frommockidtozqso[i] for i in mockiddla ]
			f_radla = sp.append(f_radla,radla)
			f_decdla = sp.append(f_decdla,decdla)
			f_zqsodla = sp.append(f_zqsodla,zqsodla)
			f_mockiddla = sp.append(f_mockiddla,mockiddla)
			f_zdla = sp.append(f_zdla,h['DLA']['Z_DLA'][:])
			f_nhidla = sp.append(f_nhidla,h['DLA']['N_HI_DLA'][:])

			h.close()
			### Sort
			w = sp.argsort(f_mockiddla)
			for el in [f_radla,f_decdla,f_mockiddla,f_zdla,f_mockiddla,f_mockiddla,f_mockiddla,f_nhidla,f_zqsodla]:
				el = el[w]

			### Save
		out = fitsio.FITS(OUTPUT,'rw',clobber=True)
		cols=[f_radla,f_decdla,f_mockiddla,f_zdla,f_mockiddla,f_mockiddla,f_mockiddla,f_nhidla,f_zqsodla]
		names = ['RA','DEC','DLAID','ZDLA','PLATE','MJD','FIBERID','NHI','ZQSO']
		out.write(cols,names=names)
		out.close()
		
		print('Done creating DLA catalog')
	return 

def create_arg_parser():
	""""Creates and returns the ArgumentParser object."""

	parser = argparse.ArgumentParser(description='Generates the DLA catalog from the quickquasars simulations output or from transmission files')
	parser.add_argument('--indir',default='./',
                    help='Path to the input directory.')
	parser.add_argument('--outdir',default='./DLA.fits.gz',
                    help='Output file')
	parser.add_argument('--transmission',action='store_true', help='use this option if the catalog will be created from the transmision file')
	

	return parser




if __name__ == "__main__":
	arg_parser = create_arg_parser()
	parsed_args = arg_parser.parse_args(sys.argv[1:])

	ini,out=parsed_args.indir,parsed_args.outdir   
	if parsed_args.transmission is True:     
		DLAS_=getDLAs_trans(ini,out)
	else: 
		DLAS_ = getDLAs_quick(ini,out)

	

