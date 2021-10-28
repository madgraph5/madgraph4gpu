part_n = 4
mome_n = 3
evnt_n = 2

# calculate aos position
evnt_i = 1
part_i = 3
mome_i = 2

aos_pos = (evnt_i * part_n * mome_n) + \
          (part_i * mome_n) + mome_i
print(aos_pos)

# calculate aosoa position
page_i = 0
evnt_i = 1
part_i = 3
mome_i = 2

aosoa_pos = (page_i * part_n * mome_n * evnt_n) + \
            (part_i * mome_n * evnt_n) + \
            (mome_i * evnt_n) + evnt_i
print(aosoa_pos)

# calculate aos elements from position
pos = 23

(evnt_i, rest) = divmod(pos, (part_n * mome_n))
(part_i, mome_i) = divmod(rest, mome_n)

print(evnt_i, part_i, mome_i)

# calculate aosoa elements from position
(page_i, rest1) = divmod(pos, (part_n * mome_n * evnt_n))
(part_i, rest2) = divmod(rest1, (mome_n * evnt_n))
(mome_i, evnt_i) = divmod(rest2, evnt_n)

print(page_i, part_i, mome_i, evnt_i)
