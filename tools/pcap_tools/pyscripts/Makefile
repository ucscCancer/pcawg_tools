CGHUB_SERVER=https://stage.cghub.ucsc.edu
TEST_UUID=40d87215-c864-49b3-8148-d5304a1a7f89
TEST_FILENAME=PAWG.pancan_bam_header.txt.bam
TEST_CHECKSUM=f6df4bbdac6add9ee4979953d569bbd0
TAR_DIR=pancan_pawg_metadata_generator

test:
	./create_pawg_metadata.py -u $(TEST_UUID) -f $(TEST_FILENAME) -c $(TEST_CHECKSUM) -p `pwd`
	./cgsubmit --validate-only -s $(CGHUB_SERVER) -u $(TEST_UUID)

tar:
	mkdir $(TAR_DIR)
	rsync -av * $(TAR_DIR)/
	rm -rf $(TAR_DIR)/$(TAR_DIR)
	rsync -av ~/cgsubmit $(TAR_DIR)/
	rsync -av ~/svn_cghub/trunk/maint/CGHWSI.py $(TAR_DIR)/
	tar -cvf $(TAR_DIR).tar $(TAR_DIR)
	gzip -9 $(TAR_DIR).tar

clean:
	rm -rf $(TEST_UUID)

clean_tar:
	rm -rf $(TAR_DIR) 
	rm $(TAR_DIR).tar
