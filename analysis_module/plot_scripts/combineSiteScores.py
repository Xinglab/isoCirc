import os,sys

if __name__ == '__main__':
    if len(sys.argv) != 11:
        print('{} xxx.score ... comb.out known/ri'.format(sys.argv[0]))
        sys.exit(1)


    if sys.argv[10] == 'known':
        out_strs = ['Known\t1\t3P\t', 'Known\t1\t5P\t',
                    'Known\t2\t3P\t', 'Known\t2\t5P\t',
                    'Novel\t1\t3P\t', 'Novel\t1\t5P\t',
                    'Novel\t2\t3P\t', 'Novel\t2\t5P\t']
        comb = sys.argv[9]
        with open(comb, 'w') as comb_fp:
            comb_fp.write('Known\tRepCnt\tType\tScores\n')
            for fn, out_str in zip(sys.argv[1:9], out_strs):
                with open(fn) as fp:
                    for line in fp:
                        comb_fp.write('{}{}\n'.format(out_str, line.rsplit()[1]))
    elif sys.argv[10] == 'ir':
        out_strs = ['Anno\tRI\t3P\t', 'Anno\tRI\t5P\t',
                    'Anno\tNonRI\t3P\t', 'Anno\tNonRI\t5P\t',
                    'All\tRI\t3P\t', 'All\tRI\t5P\t',
                    'All\tNonRI\t3P\t', 'All\tNonRI\t5P\t']
        comb = sys.argv[9]
        with open(comb, 'w') as comb_fp:
            comb_fp.write('Anno\tCate\tType\tScores\n')
            for fn, out_str in zip(sys.argv[1:9], out_strs):
                with open(fn) as fp:
                    for line in fp:
                        comb_fp.write('{}{}\n'.format(out_str, line.rsplit()[1]))
