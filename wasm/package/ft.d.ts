export interface DecodeOptions {
  /** 1-3, as jt9's -d (default 2) */
  depth?: number;
  /** Search range, Hz (default 100-3200) */
  freqLow?: number;
  freqHigh?: number;
  /** Rx and Tx audio frequencies, Hz (default 1500) */
  rxFreq?: number;
  txFreq?: number;
  /** a-priori search width, Hz (default 50) */
  apWidth?: number;
  /** Blank by default: only the CQ a-priori passes run */
  myCall?: string;
  hisCall?: string;
  hisGrid?: string;
  qsoProgress?: number;
  /** nexp_decode (its low 3 bits: the contest type) */
  contest?: number;
  /** a-priori decoding (default true) */
  apDecode?: boolean;
  apCqOnly?: boolean;
  /** Quicker subtraction: the same decodes, not WSJT-X's to the bit */
  fastSubtract?: boolean;
}

export interface Decode {
  snr: number;
  /** Seconds */
  dt: number;
  /** Audio frequency, Hz */
  freq: number;
  message: string;
  /** jt9's '?' */
  lowConfidence: boolean;
  /** a-priori type, 0 for none */
  ap: number;
  quality: number;
  sync: number;
  /** The line jt9 prints */
  line: string;
}

export interface DecodeResult {
  decodes: Decode[];
  /** <DecodeFinished>'s count */
  count: number;
  /** jt9's stdout for the period */
  text: string;
}

export interface Decoder {
  readonly mode: "ft8" | "ft4" | "ft2";
  /**
   * One period's 12 kHz audio from its start (FT8 15 s, FT4 7.5 s, FT2
   * 3.75 s; shorter is zero-filled). Float32Array is -1..1, Int16Array as is.
   * utc: HHMMSS.
   */
  decode(samples: Float32Array | Int16Array, utc: number, options?: DecodeOptions): DecodeResult;
  /** Frees the decoder's memory */
  delete(): void;
}

export interface FtModule {
  Decoder: new (mode: "ft8" | "ft4" | "ft2") => Decoder;
}

export default function createFt(moduleOverrides?: object): Promise<FtModule>;
