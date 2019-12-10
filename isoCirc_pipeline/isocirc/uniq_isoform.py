import sys
from collections import defaultdict as dd
import isocirc.utils as ut
# from .utils import *

# TODO polish uniq-isoform
# coor= [site1, site2, site3, site4 ...]
# return {read_id:isofrom_id}
def is_same_isoform(last_coor, coor, site_dis, end_dis):
	if len(last_coor) != len(coor):
		return False
	if abs(last_coor[0]-coor[0]) > end_dis: return False
	if abs(last_coor[-1]-coor[-1]) > end_dis: return False
	for i,j in zip(last_coor[1:-1], coor[1:-1]):
		if abs(i-j) > site_dis: 
			return False
	return True
	# diff = [abs(i - j) for i, j in zip(last_coor, coor)]
	# if (not diff[1:-1] or max(set(diff[1:-1])) <= site_dis) and max(set(diff[0], diff[-1])) <= end_dis:
		# return True
	# else:
		# return False

def uniq_isoform_with_unsorted_coors(coor_dict=dict(), sort_ret=False):
	ut.err_format_time('uniq_isoform_with_unsorted_coors', 'Generating isoform-wise evaluation result ... ')
	coor_to_iso_dict = dd(lambda: -1)  # name# : iso_id
	iso_to_name_dict = dd(lambda: [])  # iso_id : [name1, name2]
	iso_id = 0
	coor_list = sorted(coor_dict.items(), key=lambda k:k[1:]) if sort_ret else coor_dict.items()
	for read_id, coor in coor_list:
		if coor_to_iso_dict[tuple(coor)] != -1:
			iso_to_name_dict[coor_to_iso_dict[tuple(coor)]].append(read_id)
		else:
			coor_to_iso_dict[tuple(coor)] = iso_id
			iso_to_name_dict[iso_id] = [read_id]
			iso_id += 1
	ut.err_format_time('uniq_isoform_with_unsorted_coors', 'Generating isoform-wise evaluation result done!')
	return iso_to_name_dict


# sorted_coors = {id : [qname, rname, coors...]
def uniq_isoform_core(sorted_coor_dict=dict(), site_dis=0, end_dis=5):
	ut.err_format_time('uniq_isoform', 'Generating isoform-wise evaluation result ... ')
	name_to_iso_dict = dd(lambda: -1)  # name# : iso_id
	iso_to_name_dict = dd(lambda: [])  # iso_id : [name1, name2]

	new_isoform_id = -1
	for coor_i, coor in sorted_coor_dict.items():
		new_isoform = 1
		for i in range(coor_i - 1, -1, -1):
			if coor[0] != sorted_coor_dict[i][0] or coor[1] - sorted_coor_dict[i][1] > end_dis:
				break
			if is_same_isoform(sorted_coor_dict[i][1:], coor[1:], site_dis, end_dis):
				name_to_iso_dict[coor_i] = name_to_iso_dict[i]  # i is old
				new_isoform = 0
				break
		if new_isoform:
			new_isoform_id += 1
			name_to_iso_dict[coor_i] = new_isoform_id
		if (coor_i+1) % 100000 == 0:
			ut.err_format_time('uniq_isoform', '{} BAM records have been processed ... '.format(coor_i+1))
	for name_i, iso in name_to_iso_dict.items():
		iso_to_name_dict[iso].append(name_i)
	ut.err_format_time('uniq_isoform', 'Generating isoform-wise evaluation result done!')
	return iso_to_name_dict
