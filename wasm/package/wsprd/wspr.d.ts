export interface WsprOptions {
  /** YYMMDD and HHMM of the period's start, as wsprd takes them from a .wav's name */
  date?: string;
  time?: string;
  /** The dial, MHz (wsprd -f) */
  dialFreq?: number;
  /** 2 (WSPR-2, default) or 15 (WSPR-15) */
  mode?: 2 | 15;
  /** Fano cycles per bit (-C; default 500, as WSJT-X runs wsprd) */
  cycles?: number;
  /** OSD depth (-o; default 4, as WSJT-X; -1 for none) */
  osdDepth?: number;
  /** Deeper search (-d) */
  deep?: boolean;
  /** Quick mode (-q) */
  quick?: boolean;
  /** false: no subtraction, one pass (-s) */
  subtraction?: boolean;
  /** false: no block demodulation pass (-B) */
  blockDemod?: boolean;
  /** The stack decoder instead of Fano (-J) */
  stackDecoder?: boolean;
  /** +-150 Hz instead of +-110 (-w) */
  wide?: boolean;
  /** Fano metric bias (-z; default 0.45) */
  bias?: number;
  /** Dial error, Hz (-e) */
  dialError?: number;
  /** No drift search (-n) */
  noDrift?: boolean;
  /** false: a fresh hash table for each decode (-H) */
  hashtable?: boolean;
}

export interface WsprDecode {
  date: string;
  time: string;
  snr: number;
  /** Seconds */
  dt: number;
  /** MHz */
  freq: number;
  /** Hz per period */
  drift: number;
  /** "K1ABC FN42 37", "PJ4/K1ABC 37", "<K1ABC> FN42MI 37" */
  message: string;
  callsign: string;
  /** Empty for a type 2 message, 6 characters for type 3 */
  grid: string;
  /** dBm */
  power: number;
  sync: number;
  /** The signal's spread, Hz */
  spread: number;
  /** Decoded by the OSD */
  osd: boolean;
  /** 1-3 */
  pass: number;
  /** wsprd's stdout line, ALL_WSPR.TXT's and wspr_spots.txt's (each with its newline) */
  line: string;
  allWsprLine: string;
  spotsLine: string;
}

export interface WsprResult {
  decodes: WsprDecode[];
  /** wsprd's stdout for the period */
  text: string;
}

export interface WsprDecoder {
  /**
   * One period's 12 kHz mono audio from its start (114 s read; WSPR-15
   * 912 s; shorter is zero-filled). Int16Array, or Float32Array -1..1.
   * Throws for a period wsprd refuses.
   */
  decode(samples: Int16Array | Float32Array, options?: WsprOptions): WsprResult;
  /** hashtable.txt as wsprd writes it */
  hashtable(): string;
  /** hashtable.txt as wsprd reads it */
  loadHashtable(text: string): void;
  /** Frees the period's buffers (the hash table is kept) */
  release(): void;
  /** Frees the decoder */
  delete(): void;
}

export interface WsprModule {
  WsprDecoder: new () => WsprDecoder;
}

export default function createWspr(moduleOverrides?: object): Promise<WsprModule>;
