package com.rinke.solutions.image;

import java.awt.image.BufferedImage;

public class BufferedImageFrame {
    private final int delay;
    private final BufferedImage image;
    private final String disposal;
    private final int width, height;

    public BufferedImageFrame (BufferedImage image, int delay, String disposal, int width, int height){
        this.image = image;
        this.delay = delay;
        this.disposal = disposal;
        this.width = width;
        this.height = height;
    }

    public BufferedImageFrame (BufferedImage image){
        this.image = image;
        this.delay = -1;
        this.disposal = null;
        this.width = -1;
        this.height = -1;
    }

    public BufferedImage getImage() {
        return image;
    }

    public int getDelay() {
        return delay;
    }

    public String getDisposal() {
        return disposal;
    }

    public int getWidth() {
        return width;
    }

    public int getHeight() {
            return height;
    }
}