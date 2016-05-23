all: prep mnh ff ff-mnh

prep: 
	cd 01_prep_ideal_case && ./run_prep_ideal_case_xyz
mnh:
	cd 02_mnh             && ./run_mesonh_xyz
ff:
	cd 03_ff              && ./run_ff.sh
ff-mnh:
	cd 04_ff-mnh          && ./run_ff-mnh.sh

clean:
	cd 01_prep_ideal_case && ./clean_prep_ideal_case_xyz
	cd 02_mnh             && ./clean_mesonh_xyz
	cd 03_ff              && ./clean_ff.sh
	cd 04_ff-mnh          && ./clean_ff-mnh.sh
