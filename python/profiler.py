import pstats

if __name__ == "__main__":
    import sys
    p = pstats.Stats(sys.argv[1])
    p.strip_dirs().sort_stats('cumulative').print_stats()
