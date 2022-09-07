from scipy.io import wavfile
def main():
    samplerate, data = wavfile.read('./10000hz_sin_10min.wav')
    count = 0
    file = open("sin_wav_10min_10000hz.txt", "w")
    file.write(f"{len(data)}\n")
    # /100000 to graph easier (numbers between 0-1 mostly)
    for i in range(len(data)):
        file.write(f'{data[i][0]} {i}\n')


if __name__ == "__main__":
    main()
