#!/bin/csh -f
source $PUBSH/bin/SetUpFreeSurfer.csh 450
setenv SUBJECTS_DIR /space/md9/6/data/MMILDB/RECHARGE/BRO_LEX/MEGA/STRUCT/FSContainers/

#set recons = ( \
##'FREESURFERRECON_MT01292009_20090129.082717_1.bak' \
##'FREESURFERRECON_MP060409_20090604.150827_1' \
##'FREESURFERRECON_ds01_071019_b_20071019.103050_1' \
##'FREESURFERRECON_gp01_071102_b_20071102.165336_1XX' \
##'FREESURFERRECON_ll01_080304_20080304.142350_1' \
##'FREESURFERRECON_rj01_070920_b_20070919.150812_1' \
##'FREESURFERRECON_A01282009_20090128.172053_1' \
##'FREESURFERRECON_BM071709_20090717.154727_1' \
##'FREESURFERRECON_CC07232009_20090723.155439_1' \
##'FREESURFERRECON_CP041009_20090410.152937_1' \
##'FREESURFERRECON_GS081309_20090813.185116_1' \
##'FREESURFERRECON_HP02032009_20090203.151715_1' \
##'FREESURFERRECON_JD07212009_20090721.161508_1' \
##'FREESURFERRECON_JF07202009_20090720.171103_1' \
##'FREESURFERRECON_JF081109B_20090811.164519_1' \
##'FREESURFERRECON_JK031809TT_20090318.164639_1' \
##'FREESURFERRECON_LF081109_20090811.153654_1' \
##'FREESURFERRECON_RS07192009_20090719.151356_1' \
#'FREESURFERRECON_JM_6ys_17May09_20090517.161405_1' \
#'FREESURFERRECON_lexsem02_080311_20080311.141728_1' \
#'FREESURFERRECON_lexsem03_080409_20080409.182842_1' \
#'FREESURFERRECON_lexsem04_080416_20080416.161509_1' \
#'FREESURFERRECON_lexsem05_080506_20080506.170927_1' \
#'FREESURFERRECON_lexsem09_080527_20080527.152915_1' \
#'FREESURFERRECON_PCNDCTL810_20091221.142414_1' \
#'FREESURFERRECON_PCNDCTL811_20091221.132731_1' \
#);

set recons = ( \
'FREESURFERRECON_child_101510_P0341_20100814.120841_1' \
'FREESURFERRECON_child_021908_lexsem02_080311_scan2_20080311.141728_1' \
'FREESURFERRECON_sli_120910_SSLI_20101209.173837_1' \
);

foreach r ($recons)
cd $SUBJECTS_DIR/$r
#if (! -e $SUBJECTS_DIR/$r/bem_backup ) then
#mv -f $SUBJECTS_DIR/$r/bem $SUBJECTS_DIR/$r/bem_backup 
#mkdir bem
#endif

if (-e $SUBJECTS_DIR/$r/mri/inner_skull4.tri) then
mv -f $SUBJECTS_DIR/$r/mri/inner_skull4.tri $SUBJECTS_DIR/$r/bem
endif

cd bem
fs_surfdip $r lh
fs_surfdip $r rh
fs_surfdec $r lh
fs_surfdec $r rh
fs_viewdec lh_white_7.dec $r lh -surf inflated
fs_viewdec rh_white_7.dec $r rh -surf inflated

end