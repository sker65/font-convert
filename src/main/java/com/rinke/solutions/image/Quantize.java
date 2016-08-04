package com.rinke.solutions.image;

import static java.lang.Math.round;

import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.awt.image.DataBuffer;
import java.awt.image.DataBufferByte;
import java.awt.image.DataBufferInt;
import java.awt.image.PixelGrabber;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import javax.imageio.ImageIO;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;

public class Quantize {

	public static class RGBTriple {
		public final int[] channels;

		public RGBTriple() {
			channels = new int[3];
		}

		public RGBTriple(int color) {
			int r = (color >> 16) & 0xFF;
			int g = (color >> 8) & 0xFF;
			int b = (color >> 0) & 0xFF;
			channels = new int[] { (int) r, (int) g, (int) b };
		}

		public RGBTriple(int R, int G, int B) {
			channels = new int[] { (int) R, (int) G, (int) B };
		}
	}

	/*
	 * The authors of this work have released all rights to it and placed it in
	 * the public domain under the Creative Commons CC0 1.0 waiver
	 * (http://creativecommons.org/publicdomain/zero/1.0/).
	 * 
	 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
	 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
	 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
	 * NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
	 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
	 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
	 * USE OR OTHER DEALINGS IN THE SOFTWARE.
	 * 
	 * Retrieved from:
	 * http://en.literateprograms.org/Floyd-Steinberg_dithering_(
	 * Java)?oldid=12476
	 */
	public static class FloydSteinbergDither {
		private static int plus_truncate_uchar(int a, int b) {
			if ((a & 0xff) + b < 0)
				return 0;
			else if ((a & 0xff) + b > 255)
				return (int) 255;
			else
				return (int) (a + b);
		}

		private static int findNearestColor(RGBTriple color, RGBTriple[] palette) {
			int minDistanceSquared = 255 * 255 + 255 * 255 + 255 * 255 + 1;
			int bestIndex = 0;
			for (int i = 0; i < palette.length; i++) {
				int Rdiff = (color.channels[0] & 0xff) - (palette[i].channels[0] & 0xff);
				int Gdiff = (color.channels[1] & 0xff) - (palette[i].channels[1] & 0xff);
				int Bdiff = (color.channels[2] & 0xff) - (palette[i].channels[2] & 0xff);
				int distanceSquared = Rdiff * Rdiff + Gdiff * Gdiff + Bdiff * Bdiff;
				if (distanceSquared < minDistanceSquared) {
					minDistanceSquared = distanceSquared;
					bestIndex = i;
				}
			}
			return bestIndex;
		}

		public static int[][] floydSteinbergDither(RGBTriple[][] image, RGBTriple[] palette) {
			int[][] result = new int[image.length][image[0].length];

			for (int y = 0; y < image.length; y++) {
				for (int x = 0; x < image[y].length; x++) {
					RGBTriple currentPixel = image[y][x];
					int index = findNearestColor(currentPixel, palette);
					result[y][x] = index;

					for (int i = 0; i < 3; i++) {
						int error = (currentPixel.channels[i] & 0xff) - (palette[index].channels[i] & 0xff);
						if (x + 1 < image[0].length) {
							image[y + 0][x + 1].channels[i] = plus_truncate_uchar(image[y + 0][x + 1].channels[i], (error * 7) >> 4);
						}
						if (y + 1 < image.length) {
							if (x - 1 > 0) {
								image[y + 1][x - 1].channels[i] = plus_truncate_uchar(image[y + 1][x - 1].channels[i], (error * 3) >> 4);
							}
							image[y + 1][x + 0].channels[i] = plus_truncate_uchar(image[y + 1][x + 0].channels[i], (error * 5) >> 4);
							if (x + 1 < image[0].length) {
								image[y + 1][x + 1].channels[i] = plus_truncate_uchar(image[y + 1][x + 1].channels[i], (error * 1) >> 4);
							}
						}
					}
				}
			}
			return result;
		}

		public static void generateDither(int[] pixels, int[] p, int w, int h) {
			RGBTriple[] palette = new RGBTriple[p.length];
			for (int i = 0; i < palette.length; i++) {
				int color = p[i];
				palette[i] = new RGBTriple(color);
			}
			RGBTriple[][] image = new RGBTriple[w][h];
			for (int x = w; x-- > 0;) {
				for (int y = h; y-- > 0;) {
					int index = y * w + x;
					int color = pixels[index];
					image[x][y] = new RGBTriple(color);
				}
			}

			int[][] result = floydSteinbergDither(image, palette);
			convert(result, pixels, p, w, h);

		}

		public static void convert(int[][] result, int[] pixels, int[] p, int w, int h) {
			for (int x = w; x-- > 0;) {
				for (int y = h; y-- > 0;) {
					int index = y * w + x;
					int index2 = result[x][y];
					pixels[index] = p[index2];
				}
			}
		}
	}

	private static class PaletteColor {
		final int color;

		public PaletteColor(int color) {
			super();
			this.color = color;
		}

		@Override
		public int hashCode() {
			return color;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			PaletteColor other = (PaletteColor) obj;
			if (color != other.color)
				return false;
			return true;
		}

		public List<Integer> indices = new ArrayList<>();
	}

	public static int[] getPixels(Image image) throws IOException {
		int w = image.getWidth(null);
		int h = image.getHeight(null);
		int pix[] = new int[w * h];
		PixelGrabber grabber = new PixelGrabber(image, 0, 0, w, h, pix, 0, w);

		try {
			if (grabber.grabPixels() != true) {
				throw new IOException("Grabber returned false: " + grabber.status());
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		return pix;
	}

	/**
	 * Returns the color distance between color1 and color2
	 */
	public static float getPixelDistance(PaletteColor color1, PaletteColor color2) {
		int c1 = color1.color;
		int r1 = (c1 >> 16) & 0xFF;
		int g1 = (c1 >> 8) & 0xFF;
		int b1 = (c1 >> 0) & 0xFF;
		int c2 = color2.color;
		int r2 = (c2 >> 16) & 0xFF;
		int g2 = (c2 >> 8) & 0xFF;
		int b2 = (c2 >> 0) & 0xFF;
		return (float) getPixelDistance(r1, g1, b1, r2, g2, b2);
	}

	public static double getPixelDistance(int r1, int g1, int b1, int r2, int g2, int b2) {
		return Math.sqrt(Math.pow(r2 - r1, 2) + Math.pow(g2 - g1, 2) + Math.pow(b2 - b1, 2));
	}

	/**
	 * Fills the given fillColors palette with the nearest colors from the given
	 * colors palette until it has the given max_cols size.
	 */
	public static void fillPalette(List<PaletteColor> fillColors, List<PaletteColor> colors, int max_cols) {
		while (fillColors.size() < max_cols) {
			int index = -1;
			float minDistance = -1;
			for (int i = 0; i < colors.size(); i++) {
				PaletteColor color1 = colors.get(i);
				for (int j = 0; j < colors.size(); j++) {
					PaletteColor color2 = colors.get(j);
					if (color1 != color2) {
						float distance = getPixelDistance(color1, color2);
						if (index == -1 || distance < minDistance) {
							index = j;
							minDistance = distance;
						}
					}
				}
			}
			if( index == -1) index=0;
			PaletteColor color = colors.get(index);
			fillColors.add(color);
		}
	}

	public static void reducePaletteByAverageDistance(List<PaletteColor> colors, int max_cols, ReductionStrategy reductionStrategy) {
		while (colors.size() > max_cols) {
			int index = -1;
			float minDistance = -1;
			for (int i = 0; i < colors.size(); i++) {
				PaletteColor color1 = colors.get(i);
				float averageDistance = 0;
				int count = 0;
				for (int j = 0; j < colors.size(); j++) {
					PaletteColor color2 = colors.get(j);
					if (color1 == color2) {
						continue;
					}
					averageDistance += getPixelDistance(color1, color2);
					count++;
				}
				averageDistance /= count;
				if (minDistance == -1 || averageDistance < minDistance) {
					minDistance = averageDistance;
					index = i;
				}
			}
			PaletteColor removed = colors.remove(index);
			// find the color with the least distance:
			PaletteColor best = null;
			minDistance = -1;
			for (int i = 0; i < colors.size(); i++) {
				PaletteColor c = colors.get(i);
				float distance = getPixelDistance(c, removed);
				if (best == null || distance < minDistance) {
					best = c;
					minDistance = distance;
				}
			}
			best.indices.addAll(removed.indices);

		}
	}

	/**
	 * Reduces the given color palette until it has the given max_cols size. The
	 * colors that are closest in distance to other colors in the palette get
	 * removed first.
	 */
	public static void reducePalette(List<PaletteColor> colors, int max_cols, ReductionStrategy reductionStrategy) {
		if (reductionStrategy == ReductionStrategy.AVERAGE_DISTANCE) {
			reducePaletteByAverageDistance(colors, max_cols, reductionStrategy);
			return;
		}
		while (colors.size() > max_cols) {
			int index1 = -1;
			int index2 = -1;
			float minDistance = -1;
			for (int i = 0; i < colors.size(); i++) {
				PaletteColor color1 = colors.get(i);
				for (int j = i + 1; j < colors.size(); j++) {
					PaletteColor color2 = colors.get(j);
					if (color1 == color2) {
						continue;
					}
					float distance = getPixelDistance(color1, color2);
					if (index1 == -1 || distance < minDistance) {
						index1 = i;
						index2 = j;
						minDistance = distance;
					}
				}
			}
			PaletteColor color1 = colors.get(index1);
			PaletteColor color2 = colors.get(index2);

			switch (reductionStrategy) {
			case BETTER_CONTRAST:
				// remove the color with the lower average distance to the other
				// palette colors
				int count = 0;
				float distance1 = 0;
				float distance2 = 0;
				for (PaletteColor c : colors) {
					if (c != color1 && c != color2) {
						count++;
						distance1 += getPixelDistance(color1, c);
						distance2 += getPixelDistance(color2, c);
					}
				}
				if (count != 0 && distance1 != distance2) {
					distance1 /= (float) count;
					distance2 /= (float) count;
					if (distance1 < distance2) {
						// remove color 1;
						colors.remove(index1);
						color2.indices.addAll(color1.indices);
					} else {
						// remove color 2;
						colors.remove(index2);
						color1.indices.addAll(color2.indices);
					}
					break;
				}
				//$FALL-THROUGH$
			default:
				// remove the color with viewer mappings to the input pixels
				if (color1.indices.size() < color2.indices.size()) {
					// remove color 1;
					colors.remove(index1);
					color2.indices.addAll(color1.indices);
				} else {
					// remove color 2;
					colors.remove(index2);
					color1.indices.addAll(color2.indices);
				}
				break;
			}

		}
	}

	/**
	 * Creates an initial color palette from the given pixels and the given
	 * palette by selecting the colors with the nearest distance to the given
	 * pixels. This method also stores the indices of the corresponding pixels
	 * inside the returned PaletteColor instances.
	 */
	public static List<PaletteColor> createInitialPalette(int pixels[], int[] palette) {
		Map<Integer, Integer> used = new HashMap<>();
		ArrayList<PaletteColor> result = new ArrayList<>();

		for (int i = 0, l = pixels.length; i < l; i++) {
			double bestDistance = Double.MAX_VALUE;
			int bestIndex = -1;

			int pixel = pixels[i];
			int r1 = (pixel >> 16) & 0xFF;
			int g1 = (pixel >> 8) & 0xFF;
			int b1 = (pixel >> 0) & 0xFF;
			for (int k = 0; k < palette.length; k++) {
				int pixel2 = palette[k];
				int r2 = (pixel2 >> 16) & 0xFF;
				int g2 = (pixel2 >> 8) & 0xFF;
				int b2 = (pixel2 >> 0) & 0xFF;
				double dist = getPixelDistance(r1, g1, b1, r2, g2, b2);
				if (dist < bestDistance) {
					bestDistance = dist;
					bestIndex = k;
				}
			}

			Integer index = used.get(bestIndex);
			PaletteColor c;
			if (index == null) {
				index = result.size();
				c = new PaletteColor(palette[bestIndex]);
				result.add(c);
				used.put(bestIndex, index);
			} else {
				c = result.get(index);
			}
			c.indices.add(i);
		}
		return result;
	}

	/**
	 * Creates a simple random color palette
	 */
	public static int[] createRandomColorPalette(int num_colors) {
		Random random = new Random(101);

		int count = 0;
		int[] result = new int[num_colors];
		float add = 360f / (float) num_colors;
		for (float i = 0; i < 360f && count < num_colors; i += add) {
			float hue = i;
			float saturation = 90 + random.nextFloat() * 10;
			float brightness = 50 + random.nextFloat() * 10;
			result[count++] = Color.HSBtoRGB(hue, saturation, brightness);
		}
		return result;
	}

	public static int[] createGrayScalePalette(int count) {
		float[] grays = new float[count];
		float step = 1f / (float) count;
		grays[0] = 0;
		for (int i = 1; i < count - 1; i++) {
			grays[i] = i * step;
		}
		grays[count - 1] = 1;
		return createGrayScalePalette(grays);
	}

	/**
	 * Returns a grayscale palette based on the given shades of gray
	 */
	public static int[] createGrayScalePalette(float[] grays) {
		int[] result = new int[grays.length];
		for (int i = 0; i < result.length; i++) {
			float f = grays[i];
			result[i] = Color.HSBtoRGB(0, 0, f);
		}
		return result;
	}

	private static QuantizedImage createResultingImage(int[] pixels, List<PaletteColor> paletteColors, boolean dither, int w, int h,
				int[] colorPalette) {
		QuantizedImage res = new QuantizedImage();
		res.palette = new int[paletteColors.size()];
		res.image = new int[pixels.length];
		for (int i = 0; i < res.palette.length; i++) {
			res.palette[i] = paletteColors.get(i).color;
		}
		if (!dither) {
			for( int j = 0; j<paletteColors.size();j++) {
				PaletteColor c = paletteColors.get(j);
				for (int i : c.indices) {
					pixels[i] = c.color;
					//res.image[i] = j; 		// just store the index of color 
				}
			}
		} else {
			FloydSteinbergDither.generateDither(pixels, res.palette, w, h);
		}
		System.out.println(pixels.length+ " "+ (256*159));
		for(int i = 0; i < pixels.length; i++) {
			// for each pixel find index
			res.image[i] = findIndex(colorPalette,pixels[i]);
		}
		return res;
	}
	
	private static int findIndex(int[] colorPalette, int col) {
		for( int i = 0; i<colorPalette.length;i++) {
			if( colorPalette[i] == col) return i;
		}
		return 0;
	}

	public static class QuantizedImage {
		public int[] palette;
		public int[] image;
	}

	public static QuantizedImage quantize(int[] pixels, int widht, int heigth, int[] colorPalette, int max_cols, boolean dither, ReductionStrategy reductionStrategy) {

		// create the initial palette by finding the best match colors from the
		// given color palette
		List<PaletteColor> paletteColors = createInitialPalette(pixels, colorPalette);

		// reduce the palette size to the given number of maximum colors
		reducePalette(paletteColors, max_cols, reductionStrategy);
		assert paletteColors.size() <= max_cols;

		if (paletteColors.size() < max_cols) {
			// fill the palette with the nearest remaining colors
			List<PaletteColor> remainingColors = new ArrayList<>();
			Set<PaletteColor> used = new HashSet<>(paletteColors);
			for (int i = 0; i < colorPalette.length; i++) {
				int color = colorPalette[i];
				PaletteColor c = new PaletteColor(color);
				if (!used.contains(c)) {
					used.add(c);
					remainingColors.add(c);
				} 
			}
			fillPalette(paletteColors, remainingColors, max_cols);
		}
		assert paletteColors.size() == max_cols;

		// TODO create pixel array with palette index additionally
		
		// create the resulting image
		return createResultingImage(pixels, paletteColors, dither, widht, heigth, colorPalette);

	}

	static enum ReductionStrategy {
		ORIGINAL_COLORS, BETTER_CONTRAST, AVERAGE_DISTANCE,
	}

	static int fontColors[] = {
		0x000000, 	// black
		0x333333, 	// gray 25%
		0x666666,	// gray 50%
		0xFFFFFF,	// white
	};
	
	static int webColors[] = {
		// first vga
		0x000000, 	// black
		0x660000, 	// dark red
		0xFF0000, 	// red
		0xFF00FF, 	// purple
		0x006666, 	// teal
		0x006600, 	// green
		0x00FF00, 	// bright green
		0x00FFFF, 	// turquoise
		0x000066,	// dark blue
		0x660066, 	// violet
		0x0000FF, 	// blue
		0x333333, 	// gray 25%
		0x666666,	// gray 50%
		0x666600, 	// dark yellow
		0xFFFF00,  	// yellow
		0xFFFFFF,	// white

		// 216 web colors
		/*0x000000,*/ 0x000033, /*0x000066,*/ 0x000099, 0x0000CC, /*0x0000FF,*/ 
		0x003300, 0x003333, 0x003366, 0x003399, 0x0033CC, 0x0033FF, 
		/*0x006600,*/ 0x006633, /*0x006666,*/ 0x006699, 0x0066CC, 0x0066FF, 
		0x009900, 0x009933, 0x009966, 0x009999, 0x0099CC, 0x0099FF, 
		0x00CC00, 0x00CC33, 0x00CC66, 0x00CC99, 0x00CCCC, 0x00CCFF, 
		/*0x00FF00,*/ 0x00FF33, 0x00FF66, 0x00FF99, 0x00FFCC, /*0x00FFFF,*/ 
		0x330000, 0x330033, 0x330066, 0x330099, 0x3300CC, 0x3300FF, 
		0x333300, 0x333333, 0x333366, 0x333399, 0x3333CC, 0x3333FF, 
		0x336600, 0x336633, 0x336666, 0x336699, 0x3366CC, 0x3366FF, 
		0x339900, 0x339933, 0x339966, 0x339999, 0x3399CC, 0x3399FF, 
		0x33CC00, 0x33CC33, 0x33CC66, 0x33CC99, 0x33CCCC, 0x33CCFF, 
		0x33FF00, 0x33FF33, 0x33FF66, 0x33FF99, 0x33FFCC, 0x33FFFF, 
		/*0x660000,*/ 0x660033, /*0x660066,*/ 0x660099, 0x6600CC, 0x6600FF, 
		0x663300, 0x663333, 0x663366, 0x663399, 0x6633CC, 0x6633FF, 
		/*0x666600,*/ 0x666633, /*0x666666,*/ 0x666699, 0x6666CC, 0x6666FF, 
		0x669900, 0x669933, 0x669966, 0x669999, 0x6699CC, 0x6699FF, 
		0x66CC00, 0x66CC33, 0x66CC66, 0x66CC99, 0x66CCCC, 0x66CCFF, 
		0x66FF00, 0x66FF33, 0x66FF66, 0x66FF99, 0x66FFCC, 0x66FFFF, 
		0x990000, 0x990033, 0x990066, 0x990099, 0x9900CC, 0x9900FF, 
		0x993300, 0x993333, 0x993366, 0x993399, 0x9933CC, 0x9933FF, 
		0x996600, 0x996633, 0x996666, 0x996699, 0x9966CC, 0x9966FF, 
		0x999900, 0x999933, 0x999966, 0x999999, 0x9999CC, 0x9999FF, 
		0x99CC00, 0x99CC33, 0x99CC66, 0x99CC99, 0x99CCCC, 0x99CCFF, 
		0x99FF00, 0x99FF33, 0x99FF66, 0x99FF99, 0x99FFCC, 0x99FFFF, 
		0xCC0000, 0xCC0033, 0xCC0066, 0xCC0099, 0xCC00CC, 0xCC00FF, 
		0xCC3300, 0xCC3333, 0xCC3366, 0xCC3399, 0xCC33CC, 0xCC33FF, 
		0xCC6600, 0xCC6633, 0xCC6666, 0xCC6699, 0xCC66CC, 0xCC66FF, 
		0xCC9900, 0xCC9933, 0xCC9966, 0xCC9999, 0xCC99CC, 0xCC99FF, 
		0xCCCC00, 0xCCCC33, 0xCCCC66, 0xCCCC99, 0xCCCCCC, 0xCCCCFF, 
		0xCCFF00, 0xCCFF33, 0xCCFF66, 0xCCFF99, 0xCCFFCC, 0xCCFFFF, 
		/*0xFF0000,*/ 0xFF0033, 0xFF0066, 0xFF0099, 0xFF00CC, /*0xFF00FF,*/ 
		0xFF3300, 0xFF3333, 0xFF3366, 0xFF3399, 0xFF33CC, 0xFF33FF, 
		0xFF6600, 0xFF6633, 0xFF6666, 0xFF6699, 0xFF66CC, 0xFF66FF, 
		0xFF9900, 0xFF9933, 0xFF9966, 0xFF9999, 0xFF99CC, 0xFF99FF, 
		0xFFCC00, 0xFFCC33, 0xFFCC66, 0xFFCC99, 0xFFCCCC, 0xFFCCFF, 
		/*0xFFFF00,*/ 0xFFFF33, 0xFFFF66, 0xFFFF99, 0xFFFFCC, /*0xFFFFFF,*/ 

	};
	
	static void createWebPalette() {
		int dcol[] = { 0x00, 0x33, 0x66, 0x99, 0xCC, 0xFF };
		for( int r = 0; r<6; r++) { 
			for( int g = 0; g<6; g++) { 
				for( int b = 0; b<6; b++) {
					System.out.print(String.format("0x%06X, ", (dcol[r]<<16) | (dcol[g]<<8) | dcol[b]));
				}	
				System.out.println("");
			}	
		}	
	}
	
	static int[][] preparePalette(int[] colors) {
		int numberOfColors = colors.length;
		int rgbPalette[][] = new int[numberOfColors][5];
		for( int i = 0; i < numberOfColors; i++) {
			int col = colors[i];
			// just use 5 most significant bits for 5 target planes
			int red = (col >> 19) & 0x1F;
			int green = (col >> 11) & 0x1F;
			int blue = (col >> 3) & 0x1F;

			for( int j = 0; j < 5 ; j++ ) { // 5 target planes, build 5 color parts
				rgbPalette[i][j] =
				((blue & 1) << 2) |
				((green&1) <<1 ) |
				 (red & 1) ;
				red >>= 1; green >>= 1; blue >>= 1;  // shift right
			}
		}
		return rgbPalette;
	}
	
	private static int[] convertTo2DWithoutUsingGetRGB(BufferedImage image) {

		DataBuffer dataBuffer = image.getRaster().getDataBuffer();
		if (dataBuffer instanceof DataBufferInt) {
			int[] pixels = ((DataBufferInt) dataBuffer).getData();
			return pixels;
		}
		final int width = image.getWidth();
		final int height = image.getHeight();
		final boolean hasAlphaChannel = image.getAlphaRaster() != null;

		byte[] pixels = ((DataBufferByte) dataBuffer).getData();

		int[] result = new int[height * width];
		if (hasAlphaChannel) {
			final int pixelLength = 4;
			for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
				int argb = 0;
				argb += (((int) pixels[pixel] & 0xff) << 24); // alpha
				argb += ((int) pixels[pixel + 1] & 0xff); // blue
				argb += (((int) pixels[pixel + 2] & 0xff) << 8); // green
				argb += (((int) pixels[pixel + 3] & 0xff) << 16); // red
				result[row * width + col] = argb;
				col++;
				if (col == width) {
					col = 0;
					row++;
				}
			}
		} else {
			final int pixelLength = 3;
			for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
				int argb = 0;
				argb += -16777216; // 255 alpha
				argb += ((int) pixels[pixel] & 0xff); // blue
				argb += (((int) pixels[pixel + 1] & 0xff) << 8); // green
				argb += (((int) pixels[pixel + 2] & 0xff) << 16); // red
				result[row * width + col] = argb;
				col++;
				if (col == width) {
					col = 0;
					row++;
				}
			}
		}

		return result;
	}
	
	static String replaceExtensionTo(String newExt, String filename) {
		int p = filename.lastIndexOf(".");
		if (p != -1)
			return filename.substring(0, p) + "." + newExt;
		return filename;
	}


	
	public static void convertBmfont(String filename, String outFilename, boolean color) throws IOException {
		DataOutputStream dos = new DataOutputStream(new FileOutputStream(outFilename ));
		dos.write("BFNT".getBytes());
		dos.writeShort(1); // version
		
		FileInputStream font = new FileInputStream(filename+File.separator+"font.fnt");
		FileInputStream im = new FileInputStream(filename+File.separator+"font.png");
		
		String basename = new File(filename).getName();
		if( basename.lastIndexOf('.') != -1)basename = basename.substring(0,basename.lastIndexOf('.'));
		
		dumpFontToDataStream(font, dos, color, basename);
		//String imgName = replaceExtensionTo("png", filename);
		
		int [] colorPalette = color ? createFromInt(webColors) : createFromInt(fontColors);	
		ReductionStrategy reductionStrategy = ReductionStrategy.AVERAGE_DISTANCE;
		
		BufferedImage image = ImageIO.read(im);
		
		int[] pixels = convertTo2DWithoutUsingGetRGB(image);
//		ImageFrame original = new ImageFrame();
//		original.setImage(new File(imgName));
//		original.setTitle("Original Image");
//		original.setLocation(0, 0);

//		Image image = original.getImage();
		int width = image.getWidth(null);
		int height = image.getHeight(null);
//		int pixels[] = getPixels(image);
		
		QuantizedImage qimg = quantize(pixels, width, height, colorPalette, color?216:4, false, reductionStrategy);

		dumpImageToDataStream(dos, qimg.image, width, height, color);
		im.close();
		font.close();
		
		dos.close();
	}

	public static void main(String[] args) {
	    JFrame.setDefaultLookAndFeelDecorated(true);
	    JDialog.setDefaultLookAndFeelDecorated(true);
	    JFrame frame = new JFrame("Font Converter goDMD");
	    frame.setLayout(new FlowLayout());
	    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	    JButton button = new JButton("Select font.png B/W");
	    button.addActionListener(e->convertFont(false));
	    frame.add(button);
	    JButton button2 = new JButton("Select font.png Color");
	    button2.addActionListener(e->convertFont(true));
	    frame.add(button2);
	    frame.pack();
	    frame.setVisible(true);
	  }		

	private static void convertFont(boolean col) {
		JFileChooser fileChooser = new JFileChooser();
        int returnValue = fileChooser.showOpenDialog(null);
        if (returnValue == JFileChooser.APPROVE_OPTION) {
          File selectedFile = fileChooser.getSelectedFile();
          File dir = selectedFile.getParentFile();
		  String filename = selectedFile.getName();
		  String basename = filename.substring(0,filename.lastIndexOf('.'));
		  String outFilename = dir.getAbsolutePath()+File.separatorChar+basename+".bft";
          try {
			convertBmfont(dir.getAbsolutePath(), outFilename, col);
			} catch (IOException e) {
				e.printStackTrace();
			}
          System.out.println(selectedFile.getName());
        }
	}

	public static void convert(String args[]) throws IOException {
		
		File dir = new File("/Users/stefanri/Documents/Pinball/fonts/ConvertFonts/");
		String[] dirs = dir.list((f,name) -> new File(f,name).isDirectory() && !name.startsWith("_"));
		int k = 0;
		for(String f: dirs) {
			String conv = dir.getAbsolutePath()+File.separatorChar+f;
			//String outFilename = replaceExtensionTo("bft", conv);
			String outFilename = dir.getAbsolutePath()+File.separatorChar+"font"+k+".bft";
			System.out.println(conv);
			System.out.println(outFilename);
			convertBmfont(conv, outFilename, false);
			k++;
		}
		dir = new File("/Users/stefanri/Documents/Pinball/fonts/ConvertFonts-col/");
		dirs = dir.list((f,name) -> new File(f,name).isDirectory() && !name.startsWith("_"));
		for(String f: dirs) {
			String conv = dir.getAbsolutePath()+File.separatorChar+f;
			//String outFilename = replaceExtensionTo("bft", conv);
			String outFilename = dir.getAbsolutePath()+File.separatorChar+"font"+k+".bft";
			System.out.println(conv);
			System.out.println(outFilename);
			convertBmfont(conv, outFilename, true);
			k++;
		}
		System.exit(0);
		
		//createWebPalette();
		
		// input parameters
		String imageFileName = args[0];
		File file = new File(imageFileName);

		boolean dither = false;
		int colorPaletteSize = 216;
		int max_cols = 216;
		max_cols = Math.min(max_cols, colorPaletteSize);
		
		dumpFontAsCode("/Users/stefanri/git/bmfont/font4.fnt", "font");
			
		int[][] preparePalette = preparePalette(webColors);
		
		PrintStream out = new PrintStream("/tmp/perpal.dat");
		
		for( int i = 0; i < webColors.length; i++) {
			out.print("{ ");
			for(int j = 0; j < 5; j++) {
				out.print(String.format("\t0x%02X,",preparePalette[i][j]));
			}
			out.println("}, ");
		}
		out.close();

		// create some random color palette
		int[] colorPalette1 = createRandomColorPalette(colorPaletteSize);
		int[] colorPalette2 = createGrayScalePalette(20);

		int [] colorPalette = createFromInt(webColors);
		
		ReductionStrategy reductionStrategy = ReductionStrategy.AVERAGE_DISTANCE;

		// show the original image inside a frame
		ImageFrame original = new ImageFrame();
		original.setImage(file);
		original.setTitle("Original Image");
		original.setLocation(0, 0);

		Image image = original.getImage();
		int width = image.getWidth(null);
		int height = image.getHeight(null);
		int pixels[] = getPixels(image);
		
	/*	GifReader gifReader = new GifReader();
		List<BufferedImageFrame> gifs = gifReader.readGif(
				new FileInputStream("/Users/stefanri/Documents/Pinball/logos/cooltext192728906844032.gif"));
		int i=0;
		out = new PrintStream("/tmp/logo.txt");
		for( BufferedImageFrame gif : gifs) {
			int[] pix = convertTo2DWithoutUsingGetRGB(gif.getImage());
			QuantizedImage qimg = quantize(pix, gif.getWidth(), gif.getHeight(), 
					colorPalette, max_cols, dither, ReductionStrategy.ORIGINAL_COLORS);
			dumpImage(out, qimg.image, gif.getWidth(), gif.getHeight(), "logo"+i);
			i++;
		}
		out.close();*/
		
		QuantizedImage qimg = quantize(pixels, width, height, colorPalette, max_cols, dither, reductionStrategy);

		out = new PrintStream("/tmp/imageout.txt");
		dumpImageAsCode(out, qimg.image, width, height, "font4");
		out.close();
		
		// show the reduced image in another frame
		ImageFrame reduced = new ImageFrame();
		reduced.setImage(width, height, pixels);
		reduced.setTitle("Quantized Image (" + qimg.palette.length + " colors, dither: " + dither + ")");
		reduced.setLocation(100, 100);

	}
	
	private static void dumpFontToDataStream(InputStream is, DataOutputStream dos, boolean color, String fontname) throws IOException {
		Xml font = new Xml(is, "font");
		Xml common = font.child("common");
		int cheight = common.integer("lineHeight");
		// write 16 byte name
		byte[] name = new byte[16];
		byte[] src = fontname.getBytes();
		for( int i=0; i < 16 && i < src.length; i++) {
			name[i] = src[i];
		}
		name[15] = 0;
		dos.write(name);
		dos.writeShort(cheight);
		ArrayList<Xml> chars = font.child("chars").children("char");
		int count = 0;
		for(Xml ch: chars) if( ch.integer("id") <255 ) count++;
		dos.writeShort(count);
		for(Xml ch: chars) {
			int ascii = ch.integer("id");
			dos.write(ascii);
			if( ascii < 255 ) {
				dos.write(ch.integer("x"));
				dos.write(ch.integer("y"));
				dos.write((int)round(ch.doubleValue("width")));
				dos.write((int)round(ch.doubleValue("height")));
				dos.write(ch.integer("xoffset"));
				dos.write(ch.integer("yoffset"));
				dos.write(ch.integer("xadvance"));
			}
		}
	}

	private static void dumpFontAsCode(String filename, String name) throws IOException {
		PrintStream out = new PrintStream(filename);
		Xml font = new Xml(filename, "font");
		ArrayList<Xml> chars = font.child("chars").children("char");
		out.println("//ascii, x,y,w,h,xo,yo, xadvance");
		int count = 0;
		for(Xml ch: chars) if( ch.integer("id") <255 ) count++;
		out.println(String.format("CharMetric font4metrics[%d] = {", count));
		for(Xml ch: chars) {
			int ascii = ch.integer("id");
			if( ascii < 255 ) {
				int x = ch.integer("x");
				int xo = ch.integer("xoffset");
				int y = ch.integer("y");
				int yo =  ch.integer("yoffset");
				int w = ch.integer("width");
				int h = ch.integer("height");
				int advance = ch.integer("xadvance");
				out.println(String.format("  { '%c', %d, %d, %d, %d, %d, %d, %d }, ", ascii,x,y,w,h,xo,yo,advance));
			}
		}
		out.println("};");
		out.close();
	}
	
	private static void dumpImageToDataStream(DataOutputStream out, int[] pixels, int width, int height, boolean color ) throws IOException {
		out.writeShort(width);
		out.writeShort(height);
		out.write(color?1:0); // type
		for (int i = 0; i < pixels.length; i++) {
			out.write(pixels[i]);
		}
	}

	private static void dumpImageAsCode(PrintStream out, int[] pixels, int width, int height, String name) {
		out.println("const Raster "+name+ " =");
		out.println(" { .height = "+height+",");
		out.println("   .width = "+width+",");
		out.println("   .raster = (const uint8_t[]){ ");
		for( int y = 0; y < height; y++) {
			for(int x = 0; x< width; x++) {
				out.print(String.format("0x%02X,", pixels[y*width+x]));
			}
			out.println(" ");
		}
		out.println("} };");
	}

	private static int[] createFromInt(int[] col) {
		int[] res = new int[col.length];
		for(int i = 0; i< col.length; i++) {
			res[i] = col[i] | 0xFF000000;
		}
		return res;
	}
	
	private static void findColorMap(String name, int[] colorsNeeded, int[] colorsPal) {
		System.out.print(" { ");
		for (int i = 0; i < colorsNeeded.length; i++) {
			for (int j = 0; j < colorsPal.length; j++) {
				if( colorsNeeded[i] == colorsPal[j] ) {
					System.out.print(j);
					break;
				}
			}
			System.out.print(", ");
		}
		System.out.println(" }, //"+name);
	}
}
