import os


def suggest_threads(total=None):
    total = total or os.cpu_count() or 1
    if total > 4:
        return total - 2
    return max(1, total - 1)
